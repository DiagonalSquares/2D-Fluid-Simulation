---
name: Check local .claude directory in project root
description: Always check the local .claude directory in the project root at the start of relevant tasks
type: feedback
---

Always check the local `.claude/` directory in the project root before starting work.

**Why:** The project contains local Claude configuration (role prompt, memory) that should inform how work is approached.

**How to apply:** At the start of a conversation or before beginning a non-trivial task, check if `.claude/` exists in the project root and read the MEMORY.md index and role_prompt.md.
