# AI Assistant

DiskMINT-Nursery is the AI assistant workflow for DiskMINT. It is a companion skill for supported coding agents, including [Claude Code](https://claude.ai/code), [OpenAI Codex CLI](https://developers.openai.com/codex/cli), and compatible assistants.

**Status: available, experimental**

DiskMINT-Nursery helps users work through the full DiskMINT workflow while grounding answers in the structured reference files distributed with the DiskMINT documentation.

## What It Helps With

- Guided installation and environment verification
- Runtime assistance for model setup, parameter selection, and output interpretation
- Error diagnosis and support escalation

The skill activates when you explicitly mention DiskMINT or `import diskmint` in a conversation. This narrow activation helps avoid interference with users working on other thermochemical disk modeling codes.

## Where To Start

For full installation and usage instructions, see {doc}`../AI Features/nursery_tutorial`.

For copy-paste prompts, see {doc}`../AI Features/prompts_reference`.

For the structured files used by the assistant, see:

- {doc}`../AI Features/install_reference`
- {doc}`../AI Features/parameters_reference`
- {doc}`../AI Features/workflow_reference`
- {doc}`../AI Features/output_format_reference`
- {doc}`../AI Features/error_reference`
