# Math Summarizer

A Python tool that intelligently chunks markdown math textbooks and uses multiple LLM providers (Groq, Together AI, Anthropic, Hugging Face) to create comprehensive summaries, built with LlamaIndex.

## Features

- **Smart Chunking**: Preserves mathematical concepts, proofs, and LaTeX equations as complete units
- **Multiple LLM Providers**: Supports Groq (fast & cheap), Together AI, Anthropic Claude, and Hugging Face
- **Math-Aware Processing**: Handles LaTeX notation, theorems, and mathematical structures
- **Cost-Effective**: Groq API is ~90% cheaper than OpenAI with excellent quality
- **Coalescing Logic**: Combines individual chunk summaries into coherent, flowing output
- **CLI Interface**: Easy-to-use command-line interface with rich progress indicators

## Installation

### Prerequisites

- Python 3.9 or higher
- Poetry (for dependency management)
- API key for your chosen provider:
  - Groq API key (recommended - fast & cheap)
  - Together AI, Anthropic, or Hugging Face tokens

### Setup

1. Clone the repository:
```bash
git clone <your-repo-url>
cd math-summarizer
```

2. Install dependencies with Poetry:
```bash
poetry install
```

3. Set up environment variables:
```bash
# For Groq (recommended)
export GROQ_API_KEY="your-groq-api-key"

# Or edit .env file with your preferred provider settings
# See .env file for all configuration options
```

4. Activate the Poetry shell:
```bash
poetry shell
```

## Usage

### Basic Usage

```bash
# Uses default provider (Groq) with fastest model
math-summarizer path/to/your/textbook.md
```

### Advanced Options

```bash
# Use different providers and models
math-summarizer textbook.md \
  --provider groq \
  --model llama-3.1-70b-versatile \
  --output summary.md \
  --chunk-size 3000 \
  --temperature 0.2 \
  --verbose

# Or use Together AI
math-summarizer textbook.md \
  --provider together \
  --model mistralai/Mixtral-8x7B-Instruct-v0.1 \
  --verbose

# Or use Anthropic Claude
math-summarizer textbook.md \
  --provider anthropic \
  --model claude-3-sonnet-20240229 \
  --verbose
```

### Options

- `--output, -o`: Specify output file path (default: `input_file_summary.md`)
- `--provider`: LLM provider: groq, together, anthropic, huggingface (default: groq)
- `--model`: Provider-specific model name
- `--chunk-size`: Maximum chunk size in tokens (default: 4000)
- `--chunk-method`: Chunking method: tokens or lines (default: tokens)
- `--temperature`: Temperature for generation (default: 0.7)
- `--verbose, -v`: Enable verbose logging

## Project Structure

```
src/math_summarizer/
├── chunking/           # Smart markdown chunking
├── summarization/      # GPT-based summarization
├── coalescing/         # Summary combination logic
├── utils/              # Utility functions
└── main.py            # CLI entry point
```

## Configuration

The tool uses environment variables for configuration. See `.env` file for all available options:

- **Provider Settings**: Choose between Groq, Together AI, Anthropic, or Hugging Face
- **API Keys**: Set the appropriate API key for your chosen provider
- **Model Settings**: Configure provider-specific models
- **Chunking Settings**: Chunk sizes, overlap settings
- **Generation Settings**: Temperature, max tokens, sampling parameters
- **Output Settings**: Format preferences, math notation handling

### Cost Comparison (per 1M tokens)

| Provider | Input Cost | Speed | Quality |
|----------|------------|-------|---------|
| **Groq** | $0.27 | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| Together AI | $0.60 | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| Anthropic | $3.00 | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| OpenAI | $10-30 | ⭐⭐ | ⭐⭐⭐⭐⭐ |

## How It Works

1. **Chunking**: The markdown file is intelligently split into chunks that preserve mathematical concepts and LaTeX equations
2. **Summarization**: Each chunk is processed by your chosen LLM provider to create focused summaries
3. **Coalescing**: Individual summaries are combined into a coherent final document that maintains mathematical flow

## Development

### Running Tests

```bash
poetry run pytest
```

### Code Formatting

```bash
poetry run black src/
poetry run isort src/
```

### Type Checking

```bash
poetry run mypy src/
```

## License

MIT License - see LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests and linting
5. Submit a pull request
