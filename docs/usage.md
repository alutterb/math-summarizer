# Math Summarizer Usage Guide

## Quick Start

1. **Install dependencies**:
   ```bash
   cd math-summarizer
   poetry install
   ```

2. **Set up environment** (choose your preferred API provider):
   ```bash
   # For Groq API (Recommended: Fast & Cheap)
   export GROQ_API_KEY="your-groq-api-key"
   
   # Or for Together AI
   export TOGETHER_API_KEY="your-together-key"
   export LLM_PROVIDER="together"
   
   # Or for Anthropic Claude
   export ANTHROPIC_API_KEY="your-anthropic-key"
   export LLM_PROVIDER="anthropic"
   
   # Or edit the .env file with your preferred settings
   ```

3. **Run the summarizer**:
   ```bash
   poetry run math-summarizer examples/sample_input.md
   
   # Or specify provider explicitly
   poetry run math-summarizer examples/sample_input.md --provider groq
   ```

## Detailed Usage

### Command Line Options

```bash
math-summarizer [OPTIONS] INPUT_FILE
```

**Arguments:**
- `INPUT_FILE`: Path to the markdown file to summarize (required)

**Options:**
- `--output, -o PATH`: Output file path (default: `input_file_summary.md`)
- `--chunk-size INTEGER`: Maximum chunk size in tokens (default: 4000)
- `--chunk-method CHOICE`: Chunking method: tokens or lines (default: tokens)
- `--lines-per-chunk INTEGER`: Lines per chunk when using line-based chunking (default: 100)
- `--provider CHOICE`: LLM provider to use: groq, together, anthropic, huggingface (default: groq)
- `--model TEXT`: Model name to use (provider-specific)
- `--temperature FLOAT`: Temperature for generation (default: 0.7)
- `--verbose, -v`: Enable verbose logging
- `--mode CHOICE`: Execution mode: api (remote APIs) or local (default: api)
- `--help`: Show help message

### Examples

**Basic usage:**
```bash
math-summarizer textbook.md
```

**Custom output file:**
```bash
math-summarizer textbook.md --output my_summary.md
```

**Use different provider and model:**
```bash
# Use Groq with Llama 3.1 70B
math-summarizer textbook.md \
  --provider groq \
  --model llama-3.1-70b-versatile \
  --chunk-size 3000 \
  --temperature 0.2 \
  --verbose

# Use Together AI with Mixtral
math-summarizer textbook.md \
  --provider together \
  --model mistralai/Mixtral-8x7B-Instruct-v0.1 \
  --verbose

# Use Anthropic Claude
math-summarizer textbook.md \
  --provider anthropic \
  --model claude-3-sonnet-20240229 \
  --verbose
```

## Configuration

### Environment Variables

Create a `.env` file in the project root with the following variables:

```bash
# =============================================================================
# LLM Provider Configuration (CHOOSE ONE)
# =============================================================================

# Default provider - Groq (Recommended: Fast & Cheap)
LLM_PROVIDER=groq
GROQ_API_KEY=your-groq-api-key-here
GROQ_MODEL=llama-3.1-8b-instant

# Alternative: Together AI
# LLM_PROVIDER=together
# TOGETHER_API_KEY=your-together-api-key-here
# TOGETHER_MODEL=meta-llama/Llama-3-8b-chat-hf

# Alternative: Anthropic Claude
# LLM_PROVIDER=anthropic
# ANTHROPIC_API_KEY=your-anthropic-api-key-here
# ANTHROPIC_MODEL=claude-3-haiku-20240307

# Alternative: Hugging Face
# LLM_PROVIDER=huggingface
# HF_TOKEN=your-hf-token-here
# HF_MODEL_NAME=facebook/bart-large-cnn

# =============================================================================
# Configuration
# =============================================================================

# Execution mode: api (remote) or local
EXECUTION_MODE=api

# API Configuration
API_TIMEOUT=120
MAX_RETRIES=3

# Chunking Configuration
MAX_CHUNK_SIZE=4000
CHUNK_OVERLAP=200
MIN_CHUNK_SIZE=100

# Generation Configuration
MAX_NEW_TOKENS=512
TEMPERATURE=0.7
DO_SAMPLE=true
TOP_P=0.9
TOP_K=50

# Output Configuration
OUTPUT_FORMAT=markdown
PRESERVE_MATH_NOTATION=true

# Logging
LOG_LEVEL=INFO
LOG_FILE=logs/math_summarizer.log
```

### Input File Format

The tool expects markdown files with mathematical content. It works best with:

- **LaTeX notation**: Both inline (`$x + y$`) and display (`$$equation$$`) math
- **Structured content**: Use markdown headers (`#`, `##`, `###`) for sections
- **Mathematical structures**: Theorems, definitions, proofs, lemmas, etc.

### Example Input Structure

```markdown
# Chapter 1: Linear Algebra

## 1.1 Vector Spaces

### Definition 1.1.1
A **vector space** $V$ over a field $F$ is...

### Theorem 1.1.2
Let $V$ be a vector space. Then $\mathbf{0} \cdot \mathbf{v} = \mathbf{0}$.

**Proof**: We proceed as follows...
```

## Output Format

The tool generates a structured summary that:

1. **Preserves mathematical notation** exactly as written
2. **Maintains logical flow** between concepts
3. **Highlights key results** (theorems, definitions, etc.)
4. **Provides section organization** with clear headers
5. **Includes metadata** about the summarization process

### Sample Output Structure

```markdown
# Mathematics Textbook Summary

## Summary Information
- **Source:** textbook.md
- **Generated:** 2024-01-15T10:30:00
- **Model:** gpt-4-turbo-preview
- **Chunks processed:** 15

---

## Chapter 1: Linear Algebra

### Key Concepts
- Vector spaces and their properties
- Linear independence and basis
- Dimension theory

### Main Results
**Theorem 1.1.2**: For any vector space $V$...

[Detailed mathematical content preserving all notation]
```

## Best Practices

### Input Preparation

1. **Use clear section headers** to help the chunker organize content
2. **Keep mathematical notation consistent** throughout the document
3. **Include complete proofs** rather than partial sketches
4. **Use standard mathematical terminology** and notation

### Optimization Tips

1. **Adjust chunk size** based on your content:
   - Smaller chunks (2000-3000 tokens) for dense mathematical content
   - Larger chunks (4000-6000 tokens) for more descriptive text

2. **Choose appropriate provider and model**:
   - **Groq** (Recommended): Fastest, cheapest, good quality
     - `llama-3.1-8b-instant`: Fastest option
     - `llama-3.1-70b-versatile`: Higher quality
     - `mixtral-8x7b-32768`: Good for technical content
   - **Together AI**: Good balance of price and quality
     - `meta-llama/Llama-3-8b-chat-hf`: Balanced option
     - `mistralai/Mixtral-8x7B-Instruct-v0.1`: Good for math
   - **Anthropic Claude**: Highest quality, more expensive
     - `claude-3-haiku-20240307`: Fast and cheaper
     - `claude-3-sonnet-20240229`: Best balance
   - **Hugging Face**: Fallback option, may have reliability issues
     - `facebook/bart-large-cnn`: Good for summarization

3. **Set temperature wisely**:
   - Lower (0.1-0.3): More consistent, factual summaries
   - Higher (0.5-0.9): More creative, varied language

### Troubleshooting

**Large files taking too long?**
- Use Groq for fastest processing (10x faster than alternatives)
- Reduce chunk size to process smaller pieces
- Use `llama-3.1-8b-instant` for maximum speed
- Split large files into smaller sections

**Mathematical notation not preserved?**
- Check that `PRESERVE_MATH_NOTATION=true` in your config
- Ensure LaTeX notation is properly formatted in input

**Poor summary quality?**
- Try Anthropic Claude for highest quality
- Use larger models (e.g., `llama-3.1-70b-versatile`)
- Increase chunk overlap for better context
- Lower temperature for more focused summaries
- Try different chunking strategies

**API errors or costs too high?**
- Switch to Groq for reliable, cheap inference
- Groq is ~90% cheaper than OpenAI
- Together AI is also much cheaper than OpenAI
- Check API key is set correctly for your chosen provider

## Advanced Usage

### Programmatic Usage

```python
from math_summarizer.config import Config
from math_summarizer.chunking.chunker import Chunker
from math_summarizer.summarization.summarizer import Summarizer
from math_summarizer.coalescing.coalescer import Coalescer

# Set up configuration for Groq
import os
os.environ['GROQ_API_KEY'] = 'your-groq-api-key'
os.environ['LLM_PROVIDER'] = 'groq'

config = Config()

# Process document
chunker = Chunker(config)
chunks = chunker.chunk_file("textbook.md")

summarizer = Summarizer(config)
summaries = []
for chunk in chunks:
    summary = summarizer.summarize_chunk(chunk)
    summaries.append(summary)

coalescer = Coalescer(config)
final_summary = coalescer.coalesce_summaries(summaries, "textbook.md")
```

### Custom Chunking Strategies

The tool supports different chunking approaches:

- **Section-based** (default): Chunks by markdown sections
- **Token-based**: Simple token-count chunking
- **Concept-based**: Groups by mathematical concepts
- **Hybrid**: Combines section and concept-based approaches

## Support

For issues, questions, or contributions:

1. Check the logs in `logs/math_summarizer.log`
2. Run with `--verbose` for detailed debugging
3. Review the example files in `examples/`
4. Check the test files for expected behavior
