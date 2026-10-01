#!/usr/bin/env python3
"""Extract axisymmetric mesh and saved scalar fields from STAR-CCM+ .sim files."""
from __future__ import annotations
from starccm_sim_types import *

class StarSimReader:
    """Narrow reader for the STAR-CCM+ serialization in the supplied files."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.data = self.path.read_bytes()
        match = _HEADER_RE.match(self.data)
        if match is None:
            raise ValueError(f"{self.path} does not have the expected STAR binary header")
        self.state_position = int(match.group(1))
        state_text = self.data[self.state_position :].decode("utf-8", "strict")
        self.objects: dict[int, dict[str, Any]] = {}
        for line_number, line in enumerate(state_text.splitlines(), 1):
            if not line or line == "__eof__":
                continue
            obj = ast.literal_eval(_clean_literal(line))
            if not isinstance(obj, dict):
                continue
            self.objects[line_number - 3] = obj
        if 2 not in self.objects or self.objects[2].get("ClassName") != "star.common.Simulation":
            raise ValueError("could not locate the STAR Simulation object")

    def object(self, reference: int) -> dict[str, Any]:
        try:
            return self.objects[int(reference)]
        except KeyError as exc:
            raise KeyError(f"missing STAR object reference {reference}") from exc

    def objects_of_class(self, class_name: str) -> Iterable[tuple[int, dict[str, Any]]]:
        for reference, obj in self.objects.items():
            if obj.get("ClassName") == class_name:
                yield reference, obj

    def read_array(self, offset: int) -> tuple[dict[str, Any], np.ndarray]:
        newline = self.data.find(b"\n", int(offset), int(offset) + 4096)
        if newline < 0:
            raise ValueError(f"array header at byte {offset} is not newline terminated")
        header = ast.literal_eval(
            _clean_literal(self.data[int(offset) : newline].decode("latin1", "strict"))
        )
        if header.get("ClassName") != "Array":
            raise ValueError(f"byte {offset} does not point to a serialized Array")
        type_name = header.get("Type")
        dtype_map = {
            "Float8": np.dtype("<f8"),
            "Unsigned4": np.dtype("<u4"),
            "Integer4": np.dtype("<i4"),
            "Integer8": np.dtype("<i8"),
            "Character1": np.dtype("u1"),
        }
        if type_name not in dtype_map:
            raise ValueError(f"unsupported serialized array type {type_name!r}")
        count = int(header["nElements"])
        array = np.frombuffer(
            self.data,
            dtype=dtype_map[type_name],
            count=count,
            offset=newline + 1,
        ).copy()
        return header, array

    def _global_u4(self, object_reference: int) -> np.ndarray:
        obj = self.object(object_reference)
        if not obj.get("ClassName", "").startswith("GlobalArray<"):
            raise ValueError(f"object {object_reference} is not a GlobalArray")
        _, array = self.read_array(int(obj["dataKey"]))
        return array.astype(np.int64, copy=False)

    def partition_map(self, storage_manager_reference: int) -> PartitionMap:
        manager = self.object(storage_manager_reference)
        if manager.get("ClassName") != "DuplicateStorageManager":
            raise ValueError(f"object {storage_manager_reference} is not a storage manager")
        n_element = self._global_u4(int(manager["nElementArray"]))
        n_interior = self._global_u4(int(manager["nInteriorArray"]))
        if n_element.shape != n_interior.shape:
            raise ValueError("partition element/interior arrays have different lengths")
        _, flat_map = self.read_array(int(manager["mapHdKey"]))
        if len(flat_map) % 2:
            raise ValueError("ghost ownership map does not contain owner/index pairs")
        pairs = flat_map.astype(np.int64, copy=False).reshape(-1, 2)
        expected_ghosts = int(np.sum(n_element - n_interior))
        if len(pairs) != expected_ghosts:
            raise ValueError(
                f"ghost map has {len(pairs)} pairs, expected {expected_ghosts}"
            )
        ghosts: list[tuple[tuple[int, int], ...]] = []
        cursor = 0
        for n_all, n_int in zip(n_element, n_interior, strict=True):
            n_ghost = int(n_all - n_int)
            current = tuple((int(a), int(b)) for a, b in pairs[cursor : cursor + n_ghost])
            ghosts.append(current)
            cursor += n_ghost
        return PartitionMap(
            n_element=n_element.astype(np.int64),
            n_interior=n_interior.astype(np.int64),
            ghosts=tuple(ghosts),
        )

    def storage_array(self, storage_reference: int) -> np.ndarray:
        storage = self.object(storage_reference)
        class_name = str(storage.get("ClassName", ""))
        if class_name.startswith("SimpleStorage<"):
            _, flat = self.read_array(int(storage["dataKey"]))
            data_size = int(storage["dataSize"])
            if data_size == 0:
                return flat
            if len(flat) % data_size:
                raise ValueError(
                    f"storage {storage_reference} length {len(flat)} is not divisible "
                    f"by dataSize {data_size}"
                )
            components = len(flat) // data_size
            if components == 1:
                return flat
            return flat.reshape(data_size, components)
        if class_name.startswith("ListStorage<"):
            _, counts = self.read_array(int(storage["countKey"]))
            _, values = self.read_array(int(storage["listKey"]))
            if int(np.sum(counts)) != len(values):
                raise ValueError("list storage counts do not sum to the serialized list length")
            return np.array(
                np.split(values, np.cumsum(counts.astype(np.int64))[:-1]),
                dtype=object,
            )
        raise ValueError(f"unsupported storage class {class_name!r}")

    def region_reference(self, name: str) -> int:
        matches = [
            reference
            for reference, obj in self.objects_of_class("star.common.Region")
            if obj.get("PresentationName") == name
        ]
        if len(matches) != 1:
            raise ValueError(f"expected one region named {name!r}, found {len(matches)}")
        return matches[0]

    def fv_region(self, region_reference: int) -> tuple[int, dict[str, Any]]:
        matches = [
            (reference, obj)
            for reference, obj in self.objects_of_class("star.common.FvRegion")
            if int(obj.get("Region", -1)) == int(region_reference)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"expected one finite-volume region for {region_reference}, found {len(matches)}"
            )
        return matches[0]

    @staticmethod
    def _split_partitioned(values: np.ndarray, counts: np.ndarray) -> list[np.ndarray]:
        if int(np.sum(counts)) != len(values):
            raise ValueError(
                f"partition counts sum to {int(np.sum(counts))}, data length is {len(values)}"
            )
        edges = np.cumsum(counts.astype(np.int64))[:-1]
        return list(np.split(values, edges))

    def _add_face_vertices(
        self,
        cell_vertices: list[set[int]],
        face_manager_reference: int,
        cell_map: PartitionMap,
        vertex_map: PartitionMap,
        *,
        boundary: bool,
    ) -> None:
        manager = self.object(face_manager_reference)
        manager_map = manager.get("map", {})
        if "VertexList" not in manager_map or "FaceCellIndex" not in manager_map:
            raise ValueError(f"face manager {face_manager_reference} lacks connectivity")
        face_partition_map = self.partition_map(face_manager_reference)
        vertex_lists = self.storage_array(int(manager_map["VertexList"]))
        face_cells = np.asarray(self.storage_array(int(manager_map["FaceCellIndex"])), dtype=np.int64)
        if face_cells.ndim != 2 or face_cells.shape[1] != 2:
            raise ValueError("FaceCellIndex must have two components")
        vertex_chunks = self._split_partitioned(vertex_lists, face_partition_map.n_interior)
        cell_chunks = self._split_partitioned(face_cells, face_partition_map.n_interior)
        for partition, (partition_vertices, partition_cells) in enumerate(
            zip(vertex_chunks, cell_chunks, strict=True)
        ):
            for local_vertices, local_cells in zip(
                partition_vertices, partition_cells, strict=True
            ):
                global_vertices = {
                    vertex_map.global_index(partition, int(local_vertex))
                    for local_vertex in np.asarray(local_vertices, dtype=np.int64)
                }
                if boundary:
                    owner_cells = (int(local_cells[0]),)
                else:
                    owner_cells = (int(local_cells[0]), int(local_cells[1]))
                for local_cell in owner_cells:
                    global_cell = cell_map.global_index(partition, local_cell)
                    cell_vertices[global_cell].update(global_vertices)

    @staticmethod
    def _polygon_geometry(points_xr: np.ndarray) -> tuple[np.ndarray, float, float]:
        if len(points_xr) < 3:
            raise ValueError("a finite-volume cell has fewer than three unique vertices")
        centre = np.mean(points_xr, axis=0)
        angles = np.arctan2(points_xr[:, 1] - centre[1], points_xr[:, 0] - centre[0])
        ordered = points_xr[np.argsort(angles)]
        x = ordered[:, 0]
        r = ordered[:, 1]
        x_next = np.roll(x, -1)
        r_next = np.roll(r, -1)
        cross = x * r_next - x_next * r
        signed_twice_area = float(np.sum(cross))
        if abs(signed_twice_area) < 1.0e-20:
            raise ValueError("a finite-volume cell polygon has zero area")
        area = 0.5 * abs(signed_twice_area)
        cx = float(np.sum((x + x_next) * cross) / (3.0 * signed_twice_area))
        cr = float(np.sum((r + r_next) * cross) / (3.0 * signed_twice_area))
        volume = 2.0 * math.pi * area * abs(cr)
        return np.array([cx, abs(cr)], dtype=float), area, volume

    def extract_salt_mesh(
        self,
        region_name: str = "Salt",
        fields: Iterable[str] = (
            "Temperature",
            "Solidity",
            "VolumeFraction",
            "Density",
            "UserSpecifiedEnergySource",
        ),
    ) -> SaltMesh:
        region_reference = self.region_reference(region_name)
        fv_reference, fv_region = self.fv_region(region_reference)
        if bool(fv_region.get("is3D", True)):
            raise ValueError("the current extractor expects a 2-D axisymmetric mesh")

        cell_manager_reference = int(fv_region["cells"])
        face_manager_reference = int(fv_region["faces"])
        vertex_manager_reference = int(fv_region["vertices"])
        cell_manager = self.object(cell_manager_reference)
        vertex_manager = self.object(vertex_manager_reference)
        cell_map = self.partition_map(cell_manager_reference)
        vertex_map = self.partition_map(vertex_manager_reference)
        if cell_map.global_size != int(fv_region["CellCount"]):
            raise ValueError("cell partition map does not match FvRegion cell count")
        if vertex_map.global_size != int(fv_region["VertexCount"]):
            raise ValueError("vertex partition map does not match FvRegion vertex count")

        coord_reference = int(vertex_manager["map"]["Coord"])
        vertices = np.asarray(self.storage_array(coord_reference), dtype=float)
        if vertices.shape != (vertex_map.global_size, 3):
            raise ValueError(f"unexpected coordinate array shape {vertices.shape}")

        extracted_fields: dict[str, np.ndarray] = {}
        for field_name in fields:
            storage_reference = cell_manager.get("map", {}).get(field_name)
            if storage_reference is None:
                continue
            values = np.asarray(self.storage_array(int(storage_reference)), dtype=float)
            if values.shape != (cell_map.global_size,):
                raise ValueError(
                    f"field {field_name!r} has shape {values.shape}; "
                    f"expected {(cell_map.global_size,)}"
                )
            extracted_fields[field_name] = values

        cell_vertices: list[set[int]] = [set() for _ in range(cell_map.global_size)]
        self._add_face_vertices(
            cell_vertices,
            face_manager_reference,
            cell_map,
            vertex_map,
            boundary=False,
        )
        boundary_manager = self.object(int(fv_region["BoundaryMeshes"]))
        for fv_boundary_reference in boundary_manager.get("Keys", []):
            fv_boundary = self.object(int(fv_boundary_reference))
            if int(fv_boundary.get("FaceCount", 0)) <= 0:
                continue
            self._add_face_vertices(
                cell_vertices,
                int(fv_boundary["faces"]),
                cell_map,
                vertex_map,
                boundary=True,
            )

        cell_vertex_ids: list[np.ndarray] = []
        centroids = np.empty((cell_map.global_size, 2), dtype=float)
        areas = np.empty(cell_map.global_size, dtype=float)
        volumes = np.empty(cell_map.global_size, dtype=float)
        for cell_index, ids in enumerate(cell_vertices):
            id_array = np.array(sorted(ids), dtype=np.int64)
            points_xr = vertices[id_array, :2]
            centroid, area, volume = self._polygon_geometry(points_xr)
            centre = np.mean(points_xr, axis=0)
            angles = np.arctan2(points_xr[:, 1] - centre[1], points_xr[:, 0] - centre[0])
            id_array = id_array[np.argsort(angles)]
            cell_vertex_ids.append(id_array)
            centroids[cell_index] = centroid
            areas[cell_index] = area
            volumes[cell_index] = volume

        version = self.object(-2) if -2 in self.objects else next(
            obj for _, obj in self.objects_of_class("StarVersion")
        )
        simulation = self.object(2)
        return SaltMesh(
            source_path=self.path,
            simulation_name=str(simulation.get("name", self.path.stem)),
            version=version,
            region_name=region_name,
            fields=extracted_fields,
            vertices_xyz_m=vertices,
            cell_vertex_ids=tuple(cell_vertex_ids),
            cell_centroids_xr_m=centroids,
            cell_area_m2=areas,
            cell_volume_m3=volumes,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sim_files", nargs="+", type=Path)
    parser.add_argument("--region", default="Salt")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    summaries = []
    for sim_file in args.sim_files:
        mesh = StarSimReader(sim_file).extract_salt_mesh(args.region)
        summaries.append(mesh.summary())
    text = json.dumps(summaries, indent=2)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
        print(args.output)


if __name__ == "__main__":
    main()
