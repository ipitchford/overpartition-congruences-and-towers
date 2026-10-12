# Overpartition congruences and class-restricted towers

Anonymous · revised collection, 1.1.0-candidate · 12 October 2026.

Three linked, unrefereed computer-assisted proof candidates, preserving their separate statements and shared ingredients:

| Paper | Manuscript | Main contribution |
|---|---|---|
| A | [Modulo 7](papers/mod7/paper.pdf) · [source](papers/mod7/paper.tex) | Wei's congruence on three classes modulo 56, with a sufficient Sturm bound of 8064. |
| B | [Modulo 17, 19 and 23](papers/mod17-19-23/paper.pdf) · [source](papers/mod17-19-23/paper.tex) | Small-space proofs for 17/23 and an exact lift-before-projection cusp certificate for 19. |
| C | [Class-restricted towers](papers/towers/paper.pdf) · [source](papers/towers/paper.tex) | Exact eigenform propagation, specified additional families and the modulo-13 whole-class classification. |

Start with [AI_INDEX.md](AI_INDEX.md), [CLAIMS.json](CLAIMS.json), [ASSURANCE.md](ASSURANCE.md) and [reviews/revision-response.md](reviews/revision-response.md). The first-level modulo-17/23 cases of C are the same results as B, not extra discoveries.

## Reproduce

Use CPython 3.13 and `python -m pip install -r requirements.txt`, then:

```
python verify.py all
PYTHONOPTIMIZE=1 python verify.py all
python negative_controls.py
```

Replace `all` by `mod7`, `mod17-19-23` or `towers` for one paper. Explicit errors remain enabled in optimized Python. Generated logs and versioned receipts appear in `verification/`; generated coefficient data are reproducible outputs, not trusted inputs. See [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md).

Universal statements depend on the written modularity, cusp and Sturm arguments. Finite series alone do not prove them. No formal verification, external peer review or unaffiliated reproduction is claimed. See [CITATION_AUDIT.md](CITATION_AUDIT.md), [PROVENANCE.md](PROVENANCE.md) and [LICENSES.md](LICENSES.md).
