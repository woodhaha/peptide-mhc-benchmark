# Tool Strategy: 3-Tier Fallback Chain

## Design Principle

MedSci must operate regardless of which tools are available. Rather than hardcoding `paper-search` CLI commands (which may not be installed), MedSci probes available tools at runtime and selects the best available tier. The fallback chain ensures graceful degradation: if MCP servers are offline, drop to CLI tools; if CLI tools are missing, drop to universal WebSearch/WebFetch.

---

## Tier Detection (Capability Probe)

Run at skill activation to determine available tier:

### Probe 1: MCP Tools
Check if these MCP tools are registered (search tool catalog, not execute):
- `mcp__arxiv-mcp-server__search_papers`
- `mcp__pubmed-mcp-server__pubmed_search_articles`
- `mcp__zotero__zotero_search_items`

If ≥2 of 3 are available → **Tier 1 active**

### Probe 2: CLI Tools
```bash
which paper-search 2>/dev/null && which paper 2>/dev/null && echo "CLI_OK" || echo "CLI_MISSING"
```
If `CLI_OK` → **Tier 2 active** (used when Tier 1 unavailable)

### Fallback
If neither Tier 1 nor Tier 2 → **Tier 3** (always available: WebSearch + WebFetch)

---

## Tier 1: MCP Tools (Preferred)

### ArXiv MCP
| Action | Tool | Parameters |
|--------|------|-----------|
| Search papers | `mcp__arxiv-mcp-server__search_papers` | query, max_results=50 |
| Download paper | `mcp__arxiv-mcp-server__download_paper` | paper_id |
| Read paper | `mcp__arxiv-mcp-server__read_paper` | paper_id |

**Search strategy**: Query with topic + method keywords. Categories: q-bio, cs.LG, stat.ML. Date: last 3 years.

### PubMed MCP
| Action | Tool | Parameters |
|--------|------|-----------|
| Search articles | `mcp__pubmed-mcp-server__pubmed_search_articles` | query with MeSH terms |

**Search strategy**: Use MeSH terms for precision. Combine with publication type filters (Review, Clinical Trial). Filter by date.

### Zotero MCP
| Action | Tool | Parameters |
|--------|------|-----------|
| Search items | `mcp__zotero__zotero_search_items` | query |
| Get fulltext | `mcp__zotero__zotero_get_item_fulltext` | item_id |
| Search collections | `mcp__zotero__zotero_search_collections` | — |

**Use for**: Closed-access journals (paywalled), user's local library. Check Zotero before assuming a paper is inaccessible.

### Sci-Hub & Paper Download (Paper Access Tier)

**Sci-Hub Search Skill** (`sci-hub-search-skill`, local at `~/.claude/skills/sci-hub-search-skill/`):
| Action | Command |
|--------|---------|
| Search by DOI | `python sci_hub_search.py search --doi "<DOI>"` |
| Search by title | `python sci_hub_search.py search --title "<title>"` |
| Search by keyword | `python sci_hub_search.py search --keyword "<query>" --results N` |
| Download PDF | `python sci_hub_search.py download --doi "<DOI>" --output paper.pdf` |
| Get metadata | `python sci_hub_search.py metadata --doi "<DOI>"` |

**Paper Search MCP** (`paper-search-mcp`, openags): 23 sources with OA-First Fallback Chain:
1. Source-native download → 2. OpenAIRE/CORE/Europe PMC/PMC → 3. Unpaywall DOI resolution → 4. Sci-Hub (last resort, opt-in)

**OA-First Policy**: Always try open-access sources first (arXiv, bioRxiv, PubMed Central, Europe PMC). Sci-Hub is the LAST RESORT, enabled only when all OA paths fail. Users are responsible for legal compliance in their jurisdiction.

**Relevance to MedSci**: When Phase 2 (Collect+Verify) encounters paywalled papers, the OA-First fallback chain is invoked. Papers downloaded via Sci-Hub are stored with `[PAYWALL-BYPASS]` provenance tag in the literature matrix.

### ARS Deterministic Verification (v3.11)
When ARS plugin is loaded, its verification scripts use:
- Semantic Scholar API (`api.semanticscholar.org`)
- OpenAlex API (`api.openalex.org`)
- Crossref API (`api.crossref.org`)
- arXiv resolver (`arxiv.org/abs/<id>`)

These are invoked via ARS's Python scripts (`semantic_scholar_client.py`, etc.) or their API endpoints. The citation verification summary (Schema: `citation_verification_summary.schema.json`) is the gold standard.

---

## Tier 2: CLI Tools (paper-search / paper)

| Action | Command |
|--------|---------|
| Google web search | `paper-search google web "<query>"` |
| Semantic Scholar papers | `paper-search semanticscholar papers "<query>" --limit 10` |
| Semantic Scholar snippets | `paper-search semanticscholar snippets "<query>"` |
| Semantic Scholar citations | `paper-search semanticscholar citations <paper_id> --limit 10` |
| Semantic Scholar references | `paper-search semanticscholar references <paper_id> --limit 10` |
| PubMed search | `paper-search pubmed "<query>"` |
| Browse URL | `paper-search browse <url>` |
| Paper outline | `paper outline <arxiv_id>` |
| Paper skim | `paper skim <arxiv_id>` |
| Paper read section | `paper read <arxiv_id> <section>` |
| Paper BibTeX | `paper bibtex <arxiv_id>` |

**Note**: Tier 2 is the least likely to be available on this system (paper/paper-search are not installed by default). Tier 3 is the practical default for most sessions.

---

## Tier 3: Universal Fallback (Always Available)

### WebSearch
Used for broad discovery and topic mapping. Best for:
- Finding recent reviews and key papers
- Landscape mapping across multiple sub-topics
- Identifying key authors, databases, and benchmarks
- General web context (not just academic papers)

**Search pattern**: 2-3 parallel searches with complementary angles.

### WebFetch on APIs
Direct API calls for structured academic data:

| Source | URL Pattern | Best For |
|--------|------------|----------|
| Semantic Scholar | `https://api.semanticscholar.org/graph/v1/paper/search?query=<query>&fields=title,authors,year,citationCount,abstract,externalIds&limit=20` | Paper search with citation counts |
| Semantic Scholar detail | `https://api.semanticscholar.org/graph/v1/paper/<paper_id>?fields=title,authors,abstract,citationCount,references,citations` | Deep paper metadata |
| Semantic Scholar citations | `https://api.semanticscholar.org/graph/v1/paper/<paper_id>/citations?fields=title,year,citationCount&limit=20` | Forward citation graph |
| Semantic Scholar references | `https://api.semanticscholar.org/graph/v1/paper/<paper_id>/references?fields=title,year,citationCount&limit=20` | Backward citation graph |
| Crossref | `https://api.crossref.org/works?query.bibliographic=<title>&rows=5` | DOI and metadata verification |
| Crossref by DOI | `https://api.crossref.org/works/<DOI>` | Single-paper metadata lookup |
| PubMed E-utilities | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=<query>&retmax=20` | PubMed ID search |
| PubMed fetch | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=<PMID>&rettype=abstract` | Abstract fetch |
| arXiv API | `https://export.arxiv.org/api/query?search_query=<query>&max_results=20` | ArXiv paper search |
| arXiv abstract | `https://arxiv.org/abs/<arxiv_id>` | Paper abstract page |
| bioRxiv API | `https://api.biorxiv.org/details/<server>/<id>` | Preprint metadata |

### Direct Paper Access
For reading full papers:
- `https://pubmed.ncbi.nlm.nih.gov/<PMID>/` → Abstract + link to full text
- `https://arxiv.org/abs/<id>` → Full PDF for open-access papers
- `https://www.biorxiv.org/content/<doi>` → Preprint full text
- `https://pmc.ncbi.nlm.nih.gov/articles/<PMCID>/` → Open-access full text

---

## Fallback Decision Tree

```
Can we use MCP tools?
├── YES → Tier 1 (MCP)
│   ├── Use arxiv-mcp + pubmed-mcp for search
│   ├── Use zotero-mcp for paywalled papers
│   └── Use ARS verification scripts for citation checking
│
└── NO → Can we use CLI tools?
    ├── YES → Tier 2 (paper-search CLI)
    │   ├── Use paper-search for all searches
    │   └── Use paper for deep reading
    │
    └── NO → Tier 3 (WebSearch + WebFetch)
        ├── Use WebSearch for broad discovery
        ├── Use WebFetch on APIs for structured data
        └── Use WebFetch on article pages for deep reading
```

---

## Verification Mapping

How to verify citations at each tier:

| Metadata Field | Tier 1 | Tier 2 | Tier 3 |
|---------------|--------|--------|--------|
| DOI exists | Crossref MCP or ARS script | `paper-search browse api.crossref.org/works/<DOI>` | WebFetch `api.crossref.org/works/<DOI>` |
| Author list | MCP paper metadata | `paper-search semanticscholar details <id>` | WebFetch S2 API `/paper/<id>` |
| Journal/Year/Volume | PubMed MCP | `paper-search pubmed "<title>"` | WebFetch `eutils.ncbi.nlm.nih.gov` |
| Numeric claim | MCP fulltext + read | `paper read <id> results` | WebFetch paper page + scan |
| Directional claim | Same as numeric | Same as numeric | Same as numeric |

---

## Domain-Specific Search Patterns

For biomedical/immunology/epitope prediction (default domain):

### Tier 1 Search Queries
```
arxiv-mcp: "peptide MHC binding deep learning" + categories:q-bio
pubmed-mcp: "Epitopes"[MeSH] AND "Deep Learning"[MeSH]
zotero: search for IEDB, NetMHCpan, MHCflurry papers
```

### Tier 3 Search Queries
```
WebSearch: "deep learning peptide-MHC binding prediction 2024 2025"
WebSearch: "immunogenicity prediction neoantigen machine learning benchmark"
WebSearch: "protein language model epitope transformer GNN 2024"
WebFetch S2: /paper/search?query=peptide+MHC+binding+deep+learning&year=2023-2025
```

---

## Related Files

| File | Relationship |
|------|-------------|
| `WORKFLOW.md` | Phases reference tool tiers for each action |
| `CITATION_INTEGRITY.md` | Verification rules use the verification mapping above |
| `SKILL_CATALOG.md` | Maps ARS tools to their verification capabilities |
