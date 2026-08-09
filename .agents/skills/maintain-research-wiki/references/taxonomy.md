# Mutual-exclusive research taxonomy

Classify by the page's central object, not by every concept mentioned.

| Type | Directory | Test |
|---|---|---|
| source | Sources | Is this a faithful record of one source? |
| project | Projects | Is this organized around one research question? |
| system | Systems | Is the object a material, structure, defect, surface, or heterostructure? |
| phenomenon | Phenomena | Is the object an observed or predicted physical behavior? |
| physical-model | Physical-Models | Is it a theoretical representation or Hamiltonian? |
| computational-method | Computational-Methods | Is it a non-AI numerical/theoretical solution method? |
| ai-method | AI-Methods | Is it an ML model, training strategy, RAG, LLM, or agent design? |
| dataset | Datasets | Is it a reusable collection of data with provenance? |
| observable | Observables | Is it a measurable or computable quantity? |
| workflow | Workflows | Is it an ordered procedure with inputs, outputs, and acceptance criteria? |
| software | Software | Is it an executable package, service, scheduler, or infrastructure tool? |

Examples: MACE architecture is `ai-method`; the MACE CLI is `software`; a DFT-to-MACE sequence is `workflow`; its configurations are a `dataset`; using it to answer diffusion under strain is a `project`.
