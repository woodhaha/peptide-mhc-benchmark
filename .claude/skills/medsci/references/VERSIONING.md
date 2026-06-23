# MedSci Skill Version History

## v1.0.0 — Initial Merge (2026-06-20)

### Sources Merged
| Source | Version | Key Contribution |
|--------|---------|-----------------|
| deep-research | (single-file) | 5-phase web-first research: Broad Discovery → Narrow/Filter → Deep Read → Citation Graph → Synthesize |
| medical-imaging-review | v3.0.0 | 6-phase write-with-verify: Paradigm → Init → Collect+Verify → 3-Axis Outline → Per-Claim Write → Multi-Agent Review → Submission |
| ARS Plugin (academic-research-skills) | v3.11.1 | Full pipeline: 4 skills, 14 slash commands, 42 agents, deterministic verification (S2+OpenAlex+Crossref+arXiv), sprint contract system |

### Key Design Decisions

1. **Multi-file architecture (13 files)**: Router SKILL.md + 11 reference files, following the medical-imaging-review pattern proven against hallucination
2. **4-mode routing**: Autonomous Research (Mode 1) / ARS Delegation (Mode 2) / Deep Review (Mode 3) / Quick Distill (Mode 4) — routes to ARS for heavy academic pipeline work, executes directly for lighter research
3. **3-tier tool fallback strategy**: MCP (Tier 1) → CLI paper/paper-search (Tier 2) → WebSearch/WebFetch (Tier 3, always available)
4. **Biomedical/immunology default domain**: Pre-loaded with epitope prediction terminology, 3-axis taxonomy, core databases, and benchmarks. Domain template system for expansion to other fields
5. **Knowledge distillation protocol**: Observe → Extract → Merge → Apply cycle. Manual trigger via "distill from [skill]". Confidence-gated merging
6. **11 hallucination patterns**: 9 from MIR (adapted for biomedical domain) + 2 novel (database version drift, domain-spanning hyperbole), crossed with ARS 7-mode failure checklist
7. **Routing over replication**: ARS has 898 files — medsci routes to it, doesn't duplicate it

### Failure Modes Addressed from Predecessors

| Predecessor Failure | MedSci Fix |
|---------------------|------------|
| deep-research hardcodes paper/paper-search CLI (unavailable on this system) | TOOL_STRATEGY.md 3-tier fallback; Tier 3 WebSearch/WebFetch always works |
| medical-imaging-review v2 shipped 17 placeholder DOIs + 30-40 citation drift | CITATION_INTEGRITY.md 5-rule protocol, adapted for biomedical |
| MIR is medical-imaging specific, hard to adapt | DOMAINS.md with default biomedical domain + template system |
| ARS is complex to navigate (14 commands, 42 agents) | SKILL_CATALOG.md quick reference + Mode 2 auto-routing |
| No mechanism to absorb knowledge from other skills | DISTILLATION_PROTOCOL.md — novel capability |

---

## Iteration Triggers

| Trigger | Action |
|---------|--------|
| A draft needs factual reset after medsci production | Improve verification rules in CITATION_INTEGRITY.md |
| A new domain is adapted 3+ times | Promote from DOMAINS.md Section C to Section A |
| A tool tier changes (new MCP, CLI deprecation) | Update TOOL_STRATEGY.md |
| A new hallucination pattern discovered | Add to HALLUCINATION_PATTERNS.md, bump version |
| A skill is distilled | Log in DISTILLATION_LOG.md, update DOMAINS.md Section C |
| ARS plugin updates to new major version | Update SKILL_CATALOG.md, verify routing compatibility |
| Mode routing produces wrong mode >10% of the time | Refine trigger keywords in SKILL.md |

---

## Version Bump Checklist

- [ ] Update version number in this file
- [ ] Document what changed and why
- [ ] Update SKILL.md frontmatter version
- [ ] Update prior_version_lessons in SKILL.md metadata
- [ ] If hallucination patterns changed: cross-check with QUALITY_CHECKLIST.md
- [ ] If domain defaults changed: update DOMAINS.md Section A
- [ ] If tool strategy changed: test all 3 tiers
