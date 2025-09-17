# LLM API Setup Guide

This guide shows you how to set up different API providers for the Math Summarizer.

## 🏆 Recommended: Groq API (Fastest & Cheapest)

### 1. Get Groq API Key
- Visit: https://console.groq.com/keys
- Sign up for a free account
- Create an API key

### 2. Set Environment Variable
```bash
export GROQ_API_KEY="your-groq-api-key-here"
```

### 3. Run with Groq
```bash
# Install dependencies
poetry install

# Run with Groq (default provider)
python -m math_summarizer.main inputs/ch1.md

# Or explicitly specify Groq
python -m math_summarizer.main inputs/ch1.md --provider groq
```

**Available Groq Models:**
- `llama-3.1-8b-instant` (default, fastest)
- `llama-3.1-70b-versatile` (more capable)
- `mixtral-8x7b-32768` (good for math)

## Alternative Providers

### Together AI (Good Balance)

```bash
# Get API key from: https://api.together.xyz/settings/api-keys
export TOGETHER_API_KEY="your-together-key"

# Run with Together AI
python -m math_summarizer.main inputs/ch1.md --provider together
```

### Anthropic Claude (Best Quality)

```bash
# Get API key from: https://console.anthropic.com/
export ANTHROPIC_API_KEY="your-anthropic-key"

# Run with Claude
python -m math_summarizer.main inputs/ch1.md --provider anthropic
```

### Hugging Face (Original, but problematic)

```bash
# Get token from: https://huggingface.co/settings/tokens
export HF_TOKEN="your-hf-token"

# Run with HuggingFace
python -m math_summarizer.main inputs/ch1.md --provider huggingface --model facebook/bart-large-cnn
```

## Cost Comparison (per 1M tokens)

| Provider | Input Cost | Output Cost | Speed |
|----------|------------|-------------|-------|
| **Groq** | $0.27 | $0.27 | ⭐⭐⭐⭐⭐ |
| Together AI | $0.60 | $0.60 | ⭐⭐⭐⭐ |
| Anthropic | $3.00 | $15.00 | ⭐⭐⭐ |
| OpenAI GPT-4 | $10.00 | $30.00 | ⭐⭐ |

## Quick Start Commands

```bash
# 1. Install dependencies
poetry install

# 2. Set your API key (choose one)
export GROQ_API_KEY="your-key"        # Recommended
export TOGETHER_API_KEY="your-key"    # Alternative
export ANTHROPIC_API_KEY="your-key"   # High quality

# 3. Run the summarizer
python -m math_summarizer.main inputs/ch1.md --verbose
```

## Troubleshooting

### If you get import errors:
```bash
poetry install  # Install all dependencies
```

### If you get API key errors:
- Make sure you've set the correct environment variable
- Check that your API key is valid
- Verify you have credits/quota remaining

### If you want to use a different model:
```bash
# Groq models
python -m math_summarizer.main inputs/ch1.md --provider groq --model llama-3.1-70b-versatile

# Together AI models
python -m math_summarizer.main inputs/ch1.md --provider together --model meta-llama/Llama-3-70b-chat-hf
```
