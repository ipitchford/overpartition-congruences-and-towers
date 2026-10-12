# AI index — Overpartition congruences and class-restricted towers

## Identity and status
Anonymous, 1.1.0-candidate, DOI 10.5281/zenodo.23319157. Unrefereed computer-assisted proof collection. Start at [README.md](README.md), [CLAIMS.json](CLAIMS.json), [STATUS.md](STATUS.md). A published DOI is availability, not correctness.

## Claims and proof map
- A7: [modulo7 source](papers/mod7/paper.tex), label thm; [PDF](papers/mod7/paper.pdf). Dissection -> residue-preserving denominator -> three eta quotients -> sufficient twist level3136 -> Sturm8064.
- B1723 and B19: [middle paper](papers/mod17-19-23/paper.tex), labels thm:1723/thm:19, with Lemma3, Proposition4, Lemma5 and cusp appendix; [PDF](papers/mod17-19-23/paper.pdf). B19 needs all 240 cusp obligations before its Hecke relation.
- C-TOWERS and C-CLASS13: [tower paper](papers/towers/paper.tex), thm:main, lem:propagation, prop:class13; [PDF](papers/towers/paper.pdf). Odd-class hypotheses are not optional. Base cases overlap B.

## Evidence and commands
Install [requirements.txt](requirements.txt), read [ENVIRONMENT.txt](ENVIRONMENT.txt), run [verify.py](verify.py) with all or a paper name. [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md) gives exact counts, bounds and excluded historical scans. [negative_controls.py](negative_controls.py) requires corrupted expectations to fail in normal and optimized modes. A zero exit is only meaningful for these declared obligations.

## Trust boundaries and limitations
[ASSURANCE.md](ASSURANCE.md) separates finite computation, written modularity, formal verification, unaffiliated reproduction and review. The exact cusp engine is not independently reimplemented. The discovery scan is filtered, not an exhaustive classification. No direct applied benefit, acceleration or four-star outcome is established.

## Prior art, attribution and objections
[CITATION_AUDIT.md](CITATION_AUDIT.md) corrects the earlier square-prime comparison, credits Wei and the known eta eigenform, and maps overlap. [reviews/revision-response.md](reviews/revision-response.md) disposes of every required review item and identifies tractable extensions. Bounded searches do not prove priority.

## Integrity, provenance and licence
Use MANIFEST.json to verify shipped file hashes. [PROVENANCE.md](PROVENANCE.md), [LICENSES.md](LICENSES.md), [LICENSE](LICENSE), [LICENSE-CODE](LICENSE-CODE) distinguish original content from dependencies. Generated logs are reproducible outputs and excluded from the source manifest; private reviews and credentials are not shipped.

## Open directions
Independently implement the symbolic cusp expansion engine; establish general existence conditions for useful lift-before-projection eta factors; classify additional prime/class families with full hypotheses and prior-art checks. Neither a few numerical examples nor model agreement closes these directions.
