"""Laboratorio 03a: ensamblaje de fragmentos de ADN.

Implementa:
1. complemento reverso;
2. solapamiento sufijo-prefijo;
3. búsqueda del camino hamiltoniano de mayor solapamiento;
4. construcción de un subgrafo acíclico con un umbral de linkage.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from PIL import Image, ImageDraw, ImageFont


FRAGMENTS: Dict[str, str] = {
    "f1": "ATCCGTTGAAGCCGCGGGC",
    "f2": "TTAACTCGAGG",
    "f3": "TTAAGTACTGCCCG",
    "f4": "ATCTGTGTCGGG",
    "f5": "CGACTCCCGACACA",
    "f6": "CACAGATCCGTTGAAGCCGCGGG",
    "f7": "CTCGAGTTAAGTA",
    "f8": "CGCGGGCAGTACTT",
}


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return sequence.translate(table)[::-1]


def orientations(fragments: Dict[str, str]) -> List[Tuple[str, str]]:
    """Returns (name+orientation, sequence) for both strands."""
    result: List[Tuple[str, str]] = []
    for name, sequence in fragments.items():
        result.extend(((f"{name}+", sequence),
                       (f"{name}-", reverse_complement(sequence))))
    return result


def overlap(left: str, right: str) -> int:
    """Largest k for which suffix(left, k) == prefix(right, k)."""
    limit = min(len(left), len(right))
    for k in range(limit, 0, -1):
        if left[-k:] == right[:k]:
            return k
    return 0


@dataclass(frozen=True)
class State:
    score: int
    path: Tuple[int, ...]


def assemble(fragments: Dict[str, str], target_length: int = 55):
    """Finds the max-overlap Hamiltonian path.

    Each original fragment is used exactly once, while either orientation
    may be selected. Ties are resolved by distance to target_length.
    """
    names = list(fragments)
    oriented = orientations(fragments)
    # (mask, final oriented fragment) -> best overlap score and path
    dp: Dict[Tuple[int, int], State] = {}
    for j, (label, _) in enumerate(oriented):
        original = names.index(label[:-1])
        dp[(1 << original, j)] = State(0, (j,))

    for size in range(1, len(names)):
        current = list(dp.items())
        for (mask, last), state in current:
            if mask.bit_count() != size:
                continue
            left = oriented[last][1]
            for nxt, (label, right) in enumerate(oriented):
                original = names.index(label[:-1])
                if mask & (1 << original):
                    continue
                new_key = (mask | (1 << original), nxt)
                new_state = State(state.score + overlap(left, right),
                                  state.path + (nxt,))
                previous = dp.get(new_key)
                if previous is None or new_state.score > previous.score:
                    dp[new_key] = new_state

    full_mask = (1 << len(names)) - 1
    candidates = []
    for (mask, _), state in dp.items():
        if mask != full_mask:
            continue
        path = state.path
        sequence = oriented[path[0]][1]
        links = []
        for left_index, right_index in zip(path, path[1:]):
            link = overlap(oriented[left_index][1], oriented[right_index][1])
            links.append(link)
            sequence += oriented[right_index][1][link:]
        candidates.append((abs(len(sequence) - target_length),
                           -state.score, sequence, path, links))

    _, _, consensus, path, links = min(candidates)
    return oriented, consensus, path, links


def dag_edges(oriented, path, threshold: int):
    """Edges of the acyclic graph induced by the selected path order."""
    order = {index: position for position, index in enumerate(path)}
    edges = []
    selected = set(path)
    for left_index, (_, left) in enumerate(oriented):
        for right_index, (_, right) in enumerate(oriented):
            if (left_index not in selected or right_index not in selected
                    or left_index == right_index
                    or order[left_index] >= order[right_index]):
                continue
            link = overlap(left, right)
            if link >= threshold:
                edges.append((oriented[left_index][0], oriented[right_index][0], link))
    return edges


def save_graph_image(oriented, path, threshold: int, output_file: str = "grafo_ensamblaje.png"):
    """Saves a publication-ready PNG of the acyclic overlap graph."""
    edges = dag_edges(oriented, path, threshold)
    labels = [oriented[index][0] for index in path]
    width, height = 2200, 650
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    try:
        title_font = ImageFont.truetype("arial.ttf", 30)
        label_font = ImageFont.truetype("arial.ttf", 23)
        value_font = ImageFont.truetype("arial.ttf", 20)
    except OSError:
        title_font = label_font = value_font = ImageFont.load_default()

    margin, y, step = 125, 360, 275
    position = {label: (margin + i * step, y) for i, label in enumerate(labels)}
    selected = {(labels[i], labels[i + 1]) for i in range(len(labels) - 1)}
    for left, right, link in edges:
        x1, y1 = position[left][0], position[left][1]
        x2, y2 = position[right][0], position[right][1]
        is_selected = (left, right) in selected
        color = "#1f4e79" if is_selected else "#9aa7b2"
        line_width = 6 if is_selected else 3
        if x2 - x1 > step + 5:
            arc_y = y - 85 - 35 * min(x2 - x1, 5) / step
            draw.arc((x1, arc_y, x2, y + 25), 180, 360, fill=color, width=line_width)
            end_y = y - 3
            draw.polygon([(x2, end_y), (x2 - 18, end_y - 13), (x2 - 13, end_y + 8)], fill=color)
            tx, ty = (x1 + x2) // 2, arc_y + 10
        else:
            draw.line((x1 + 45, y, x2 - 45, y), fill=color, width=line_width)
            draw.polygon([(x2 - 35, y), (x2 - 55, y - 13), (x2 - 55, y + 13)], fill=color)
            tx, ty = (x1 + x2) // 2, y - 45
        text_color = "#b22222" if is_selected else "#58636e"
        draw.text((tx, ty), str(link), font=value_font, fill=text_color, anchor="mm")

    draw.text((width // 2, 55), f"Subgrafo acíclico de solapamientos (linkage t = {threshold})",
              font=title_font, fill="#1f4e79", anchor="mm")
    for label, (x, node_y) in position.items():
        draw.ellipse((x - 45, node_y - 45, x + 45, node_y + 45), fill="#eaf2f8", outline="#1f4e79", width=5)
        draw.text((x, node_y), label, font=label_font, fill="#1f1f1f", anchor="mm")
    draw.text((width // 2, 585), "Número sobre cada arista = longitud del solapamiento",
              font=value_font, fill="#444444", anchor="mm")
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    image.save(output_file)


if __name__ == "__main__":
    oriented, consensus, path, links = assemble(FRAGMENTS, target_length=55)
    selected = [oriented[index][0] for index in path]
    print("Camino:", " -> ".join(selected))
    print("Linkages:", links)
    print("Longitud:", len(consensus))
    print("Consenso:", consensus)
    threshold = min(links)
    print(f"Umbral linkage: {threshold}")
    save_graph_image(oriented, path, threshold, "outputs/grafo_ensamblaje.png")
    print("Imagen del grafo: outputs/grafo_ensamblaje.png")
    print("Aristas del DAG:")
    for left, right, link in dag_edges(oriented, path, threshold):
        print(f"  {left} -> {right} ({link})")
