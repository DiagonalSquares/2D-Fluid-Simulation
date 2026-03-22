---
name: feedback_local_memory
description: User wants project knowledge stored in .claude/ within the project, not in global persistent memory
type: feedback
---

Store project knowledge bases and memory files in the project's local `.claude/` directory, NOT in the global `~/.claude/projects/` persistent memory path.

**Why:** User explicitly corrected this — they want project knowledge colocated with the project.
**How to apply:** When creating memory or knowledge files for this project, always use `.claude/memory/` relative to the project root.
