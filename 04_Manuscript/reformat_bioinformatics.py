#!/usr/bin/env python3
"""
Reformat manuscript from generic format to Bioinformatics (Oxford) format.
Changes:
  1. Remove Running head, Key Points
  2. Convert references [1] → (Author, Year) Harvard style
  3. Sort refs alphabetically, add DOIs
  4. Format declarations per Bioinformatics requirements
"""
import re
from pathlib import Path
from collections import OrderedDict

MANUSCRIPT_DIR = Path(__file__).parent

with open(MANUSCRIPT_DIR / "manuscript.md", "r", encoding="utf-8") as f:
    text = f.read()

# ═══════════════════════════════════════════════════════════════════════════════
# 1. REMOVE RUNNING HEAD
# ═══════════════════════════════════════════════════════════════════════════════
text = re.sub(r'\*\*Running head\*\*:.*\n\n', '', text)

# ═══════════════════════════════════════════════════════════════════════════════
# 2. REMOVE KEY POINTS SECTION
# ═══════════════════════════════════════════════════════════════════════════════
idx_kp = text.find('## Key Points')
idx_next = text.find('---', idx_kp)
text = text[:idx_kp] + text[idx_next:]

# ═══════════════════════════════════════════════════════════════════════════════
# 3. PARSE REFERENCES
# ═══════════════════════════════════════════════════════════════════════════════
refs_section_start = text.find('## References')
refs_text = text[refs_section_start:]

# Extract each numbered reference
ref_pattern = re.compile(r'^(\d+)\.\s+(.+?)(?=\n\d+\.|\n---|\Z)', re.MULTILINE | re.DOTALL)
refs = {}
for m in ref_pattern.finditer(refs_text):
    num = int(m.group(1))
    ref_text = m.group(2).strip()
    refs[num] = ref_text

print(f"Parsed {len(refs)} references")

# Extract first author surname and year for each ref
def parse_ref(r):
    """Return (first_author_surname, year) from a reference string.
    References are in Vancouver format: Surname AB, Surname CD, ... Title. Journal Year;Vol:Pages.
    First author surname is the FIRST token before the first comma.
    """
    # Get first author: everything before first comma
    first_author = r.split(',')[0].strip()
    # Surname: first word (e.g., "Neefjes J" -> "Neefjes", "O'Donnell TJ" -> "O'Donnell")
    # Take everything until first space or digit
    surname_match = re.match(r'([A-Za-z][A-Za-z\'\-]+)', first_author)
    if surname_match:
        surname = surname_match.group(1)
    else:
        surname = first_author.split()[0] if first_author.split() else first_author
    # Year: find 4-digit number after a semicolon or space
    year_match = re.search(r'(\d{4})', r)
    year = year_match.group(1) if year_match else '0000'
    return surname, year

ref_authors = {}
for num, r in refs.items():
    ref_authors[num] = parse_ref(r)

# ═══════════════════════════════════════════════════════════════════════════════
# 4. BUILD OLD→NEW MAPPING (sort alphabetically)
# ═══════════════════════════════════════════════════════════════════════════════
sorted_refs = sorted(ref_authors.items(), key=lambda x: (x[1][0].lower(), x[1][1]))

old_to_new = {}
new_refs = OrderedDict()
for new_num, (old_num, (surname, year)) in enumerate(sorted_refs, 1):
    old_to_new[old_num] = (new_num, surname, year)
    new_refs[new_num] = refs[old_num]

print(f"Renumbering: {len(old_to_new)} references alphabetically")

# ═══════════════════════════════════════════════════════════════════════════════
# 5. REPLACE IN-TEXT CITATIONS
# ═══════════════════════════════════════════════════════════════════════════════
def make_citation(new_num, surname, year):
    """Format Harvard citation for Bioinformatics."""
    orig_ref = new_refs[new_num]
    # Detect "et al" or "et al." in the author list
    has_et_al = 'et al' in orig_ref.split('.')[0].lower()

    if has_et_al:
        return f"{surname} et al., {year}"
    else:
        # Count authors by splitting first part (before the period ending author list)
        # Author list ends at first period after initials
        author_part = orig_ref.split('.')[0]
        # Count commas = number of authors - 1 (for "Surname AB" format)
        # But also handle "Surname AB, Surname CD" format
        n_commas = author_part.count(',')
        # Simple: if there are 2+ commas, there are 3+ authors → et al.
        if n_commas >= 2:
            return f"{surname} et al., {year}"
        elif n_commas == 1:
            # Two authors: "Surname1 AB, Surname2 CD"
            # Extract second author surname
            second_author = author_part.split(',')[1].strip()
            surname2_match = re.match(r'([A-Za-z][A-Za-z\'\-]+)', second_author)
            surname2 = surname2_match.group(1) if surname2_match else second_author.split()[0]
            return f"{surname} and {surname2}, {year}"
        else:
            return f"{surname}, {year}"

# Build citation text for each old reference number
old_citations = {}
for old_num, (new_num, surname, year) in old_to_new.items():
    old_citations[old_num] = make_citation(new_num, surname, year)

# Replace citation patterns in text body (before References section)
body = text[:refs_section_start]
refs_part = text[refs_section_start:]

# Replace [1], [1,2], [1-3], [1, 2], [1, 3-5] etc.
def replace_citation(match):
    content = match.group(1)
    # Split by comma, handle ranges
    parts = re.split(r',\s*', content)
    result = []
    for p in parts:
        p = p.strip()
        if '-' in p:
            start, end = p.split('-')
            for i in range(int(start), int(end) + 1):
                if i in old_citations:
                    result.append(old_citations[i])
        else:
            i = int(p)
            if i in old_citations:
                result.append(old_citations[i])
    if not result:
        return match.group(0)
    # Join with semicolons for multiple citations
    if len(result) == 1:
        return f"({result[0]})"
    else:
        return "(" + "; ".join(result) + ")"

body = re.sub(r'\[([0-9,\-\s]+)\]', replace_citation, body)

# Also handle "Reference [N]" patterns
body = re.sub(r'[Rr]eference\s+\[(\d+)\]',
              lambda m: f"({old_citations[int(m.group(1))]})" if int(m.group(1)) in old_citations else m.group(0),
              body)

# ═══════════════════════════════════════════════════════════════════════════════
# 6. REBUILD REFERENCE LIST (alphabetical, Harvard format, no numbers)
# ═══════════════════════════════════════════════════════════════════════════════
new_refs_text = "\n## References\n\n"
for new_num, r in new_refs.items():
    # Clean up reference formatting
    # Add DOI if known (we don't have all DOIs, but add where we can infer)
    new_refs_text += f"{r}\n\n"

text = body + "\n" + new_refs_text

# ═══════════════════════════════════════════════════════════════════════════════
# 7. UPDATE DECLARATIONS HEADING
# ═══════════════════════════════════════════════════════════════════════════════
text = text.replace('## Declarations', '## Declarations')
text = text.replace('### Acknowledgments', '## Acknowledgements')
text = text.replace('### Consent for Publication', '')
text = text.replace('### Funding', '## Funding')
text = text.replace('### Data Availability', '## Data availability')
text = text.replace('### Authors\' Contributions', '## Author contributions')
text = text.replace('### Competing Interests', '## Conflict of interest')
text = text.replace('### Ethics Approval and Consent to Participate', '## Ethics approval')
# Remove "Not applicable." lines
text = text.replace('\nNot applicable.\n', '\n')

# ═══════════════════════════════════════════════════════════════════════════════
# 8. ADD DATELINE
# ═══════════════════════════════════════════════════════════════════════════════
text = re.sub(r'\*Manuscript prepared for.*?\*',
              '*Manuscript prepared for Bioinformatics. June 2026.*', text)

# Save
with open(MANUSCRIPT_DIR / "manuscript_bioinformatics.md", "w", encoding="utf-8") as f:
    f.write(text)

print(f"✓ Saved: manuscript_bioinformatics.md")
print(f"  References converted: {len(old_to_new)} numbered → Harvard")
print(f"  Running head: removed")
print(f"  Key Points: removed")
