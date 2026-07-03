---
description: Git branching, committing, and worktree workflow.
---

# Git Workflow

The MEFS repository uses a standard git workflow centered on the `main` branch.

## 1. Branching Strategy
- All development should occur on feature branches branched from `main`.
- Merge feature branches into `main` using Pull Requests after compiling and verifying all changes.

## 2. Commit Guidelines
- Write clear, descriptive commit messages.
- Format: `<component>: <short summary of changes>`
- Example: `core: optimize FFT calculation scaling`

## 3. Worktree Usage
For maintaining multiple branches or reviewing PRs, use Git worktrees rather than cloning multiple times:
```bash
git worktree add ../mefs-review feature-branch
```
