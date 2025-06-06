import collections
import functools
import pathlib
from typing import Iterable, Literal, TYPE_CHECKING, Union, Callable

from gen._chunk._builder import Builder
from gen._chunk._complex_util import complex_key
from gen._chunk._flow import Flow, PauliMap
from gen._chunk._keyed_pauli_map import KeyedPauliMap
from gen._chunk._patch import Patch
from gen._chunk._stabilizer_code import StabilizerCode
from gen._chunk._tile import Tile

if TYPE_CHECKING:
    from gen._chunk._chunk import Chunk


class ChunkInterface:
    """Specifies a set of stabilizers and observables that a chunk can consume or prepare."""

    def __init__(
        self,
        ports: Iterable[PauliMap | KeyedPauliMap],
        *,
        discards: Iterable[PauliMap | KeyedPauliMap] = (),
    ):
        self.ports = frozenset(ports)
        self.discards = frozenset(discards)

    def with_transformed_coords(
        self, transform: Callable[[complex], complex]
    ) -> "ChunkInterface":
        return ChunkInterface(
            ports=[port.with_transformed_coords(transform) for port in self.ports],
            discards=[
                discard.with_transformed_coords(transform) for discard in self.discards
            ],
        )

    @functools.cached_property
    def used_set(self) -> frozenset[complex]:
        return frozenset(
            q for port in self.ports | self.discards for q in port.qubits.keys()
        )

    def write_svg(
        self,
        path: str | pathlib.Path,
        *,
        show_order: bool | Literal["undirected", "3couplerspecial"] = False,
        show_measure_qubits: bool = False,
        show_data_qubits: bool = True,
        system_qubits: Iterable[complex] = (),
        opacity: float = 1,
        show_coords: bool = True,
        show_obs: bool = True,
        other: Union[
            None, "StabilizerCode", "Patch", Iterable[Union["StabilizerCode", "Patch"]]
        ] = None,
        tile_color_func: Callable[["Tile"], str] | None = None,
        rows: int | None = None,
        cols: int | None = None,
        find_logical_err_max_weight: int | None = None,
    ) -> None:
        flat = [self]
        if isinstance(other, (StabilizerCode, Patch)):
            flat.append(other)
        elif other is not None:
            flat.extend(other)

        from gen._viz_patch_svg import patch_svg_viewer
        from gen._util import write_file

        viewer = patch_svg_viewer(
            patches=flat,
            show_obs=show_obs,
            show_measure_qubits=show_measure_qubits,
            show_data_qubits=show_data_qubits,
            show_order=show_order,
            find_logical_err_max_weight=find_logical_err_max_weight,
            system_qubits=system_qubits,
            opacity=opacity,
            show_coords=show_coords,
            tile_color_func=tile_color_func,
            cols=cols,
            rows=rows,
        )
        write_file(path, viewer)

    def without_discards(self) -> "ChunkInterface":
        return self.with_edits(discards=())

    def without_keyed(self) -> "ChunkInterface":
        return ChunkInterface(
            ports=[port for port in self.ports if isinstance(port, PauliMap)],
            discards=[
                discard for discard in self.discards if isinstance(discard, PauliMap)
            ],
        )

    def with_discards_as_ports(self) -> "ChunkInterface":
        return self.with_edits(discards=(), ports=self.ports | self.discards)

    def with_anonymized_keys(self) -> "ChunkInterface":
        return self.with_edits(
            ports=[
                port.pauli_string if isinstance(port, KeyedPauliMap) else port
                for port in self.ports
            ],
            discards=[
                discard.pauli_string if isinstance(discard, KeyedPauliMap) else discard
                for discard in self.discards
            ],
        )

    def __repr__(self) -> str:
        lines = ["gen.ChunkInterface("]

        lines.append(f"    ports=["),
        for port in sorted(self.ports):
            lines.append(f"        {port!r},")
        lines.append("    ],")

        if self.discards:
            lines.append(f"    discards=["),
            for discard in sorted(self.discards):
                lines.append(f"        {discard!r},")
            lines.append("    ],")

        lines.append(")")
        return "\n".join(lines)

    def __str__(self) -> str:
        lines = []
        for port in sorted(self.ports):
            lines.append(str(port))
        for discard in sorted(self.discards):
            lines.append(f"discard {discard}")
        return "\n".join(lines)

    def with_edits(
        self,
        *,
        ports: Iterable[PauliMap | KeyedPauliMap] | None = None,
        discards: Iterable[PauliMap | KeyedPauliMap] | None = None,
    ) -> "ChunkInterface":
        return ChunkInterface(
            ports=self.ports if ports is None else ports,
            discards=self.discards if discards is None else discards,
        )

    def __eq__(self, other):
        if not isinstance(other, ChunkInterface):
            return NotImplemented
        return self.ports == other.ports and self.discards == other.discards

    @functools.cached_property
    def data_set(self) -> frozenset[complex]:
        return frozenset(
            q
            for pauli_string_list in [self.ports, self.discards]
            for ps in pauli_string_list
            for q in ps.qubits
        )

    def _mpp_chunk(self, *, direction: Literal["in", "out"]) -> "Chunk":
        builder = Builder.for_qubits(self.data_set)

        layers = collections.defaultdict(list)
        qubit_used = set()
        for port in sorted(self.ports):
            layer = 0
            while any((q, layer) in qubit_used for q in port.qubits.keys()):
                layer += 1
            for q in port.qubits.keys():
                qubit_used.add((q, layer))
            layers[layer].append(port)
        for _, layer in sorted(layers.items()):
            for port in layer:
                builder.append("MPP", [port])
        flows = []
        for port in sorted(self.ports):
            flows.append(
                Flow(
                    start=port.pauli_string if direction == "in" else None,
                    end=port.pauli_string if direction == "out" else None,
                    center=min(port.qubits, key=complex_key),
                    measurement_indices=builder.lookup_recs([port]),
                    obs_key=port.key if isinstance(port, KeyedPauliMap) else None,
                )
            )

        from gen._chunk._chunk import Chunk

        return Chunk(
            circuit=builder.circuit,
            q2i=builder.q2i,
            flows=flows,
            discarded_inputs=self.discards if direction == "in" else (),
            discarded_outputs=self.discards if direction == "out" else (),
        )

    def mpp_init_chunk(self) -> "Chunk":
        return self._mpp_chunk(direction="out")

    def mpp_end_chunk(self) -> "Chunk":
        return self._mpp_chunk(direction="in")

    def to_patch(self) -> Patch:
        return Patch(
            tiles=[
                Tile(
                    bases="".join(port.qubits.values()),
                    data_qubits=port.qubits.keys(),
                    measure_qubit=min(port.qubits, key=complex_key),
                )
                for pauli_string_list in [self.ports, self.discards]
                for port in pauli_string_list
                if isinstance(port, PauliMap)
            ]
        )

    def to_code(self) -> StabilizerCode:
        return StabilizerCode(
            stabilizers=self.to_patch(),
            logicals=[
                port.pauli_string
                for pauli_string_list in [self.ports, self.discards]
                for port in pauli_string_list
                if isinstance(port, KeyedPauliMap)
            ],
        )
