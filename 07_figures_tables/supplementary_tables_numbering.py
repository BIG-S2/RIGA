"""Render Nature Genetics table numbering from stable analysis table sources.

Run from any directory with Python 3. Source CSV names remain unchanged.
This script regenerates ng_supp_tables only; manuscript citations are edited
separately when the order of first citation changes.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
ORDER = [1, 2, 8, 13, 14, 12, 10, 3, 4, 9, 11, 6, 7, 5, 15, 16, 17, 18, 19, 20]
NUMBER = {old: new for new, old in enumerate(ORDER, 1)}
TABLE_CITATION = re.compile(
    r"(Supplementary Tables?[ ~]+)"
    r"(\d+[a-z]?(?:(?:,\s*(?:and[ ~]+)?| +and[ ~]+|--)[~ ]*\d+[a-z]?)*)"
)


def renumber(text):
    text = TABLE_CITATION.sub(
        lambda m: m[1] + re.sub(r"\d+", lambda n: str(NUMBER[int(n[0])]), m[2]),
        text,
    )
    return re.sub(r"tab:s(\d+)", lambda m: "tab:s" + str(NUMBER[int(m[1])]), text)


def spell_out_test_name(text):
    """Expand reader-facing names while leaving analysis identifiers unchanged."""
    headers = {
        "SB (primary)": r"\shortstack{SKAT--burden\\omnibus test\\(primary)}",
        "$q$ (SB)": r"\shortstack{$q$ (SKAT--burden\\omnibus test)}",
        "ADSP SB $P$": r"\shortstack{ADSP SKAT--burden\\omnibus test $P$}",
        "SB $P$, conditional": r"\shortstack{SKAT--burden\\omnibus test $P$,\\conditional}",
        "SB $P$": r"\shortstack{SKAT--burden\\omnibus test $P$}",
    }
    for old, new in headers.items():
        text = text.replace(r"\textbf{" + old + "}", r"\textbf{" + new + "}")
    text = text.replace("SB (SKAT + burden omnibus)", "SKAT--burden omnibus test")
    text = text.replace("statistic SB", "SKAT--burden omnibus test")
    text = text.replace("adjusted SB", "adjusted SKAT--burden omnibus test $P$ values")
    text = text.replace("SB is the Cauchy", "The SKAT--burden omnibus test is the Cauchy")
    text = text.replace("SB is the primary", "the SKAT--burden omnibus test is the primary")
    text = text.replace("brings SB and", "brings the SKAT--burden omnibus test and")
    text = text.replace("omnibus and SB,", "omnibus and the SKAT--burden omnibus test,")
    return re.sub(r"\bSB\b", "SKAT--burden omnibus test", text)


def render():
    output = ROOT / "ng_supp_tables"
    output.mkdir(exist_ok=True)
    count = 0
    for source in sorted((ROOT / "supp_tables").glob("S*.tex")):
        match = re.fullmatch(r"S(\d+)([a-z]?)\.tex", source.name)
        if match is None or int(match[1]) not in NUMBER:
            continue
        target = output / f"S{NUMBER[int(match[1])]}{match[2]}.tex"
        text = spell_out_test_name(renumber(source.read_text()))
        text = text.replace("Proposition 1b", r"Proposition~\ref{prop:strat}(b)")
        text = text.replace("GRCh38; sequenced", "GRCh38, sequenced")
        text = text.replace("stratified design; Proposition", "stratified design, Proposition")
        text = text.replace("Analytic set -- cases", "Cases in the analytic set")
        text = text.replace("Analytic set -- controls", "Controls in the analytic set")
        target.write_text(text)
        count += 1
    mapping = "Original analysis table -> Nature Genetics table\n"
    mapping += "\n".join(f"S{old} -> S{NUMBER[old]}" for old in ORDER) + "\n"
    (output / "numbering_map.txt").write_text(mapping)
    print(f"Rendered {count} table parts in {output}")


if __name__ == "__main__":
    render()
