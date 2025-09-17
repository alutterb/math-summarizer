# Migration Guide: Multi-Provider LLM Support

This document summarizes the changes made to support multiple LLM providers (Groq, Together AI, Anthropic, Hugging Face) as alternatives to OpenAI.

## 🚀 Quick Start (New Users)

```bash
# 1. Get a Groq API key (recommended - fast & cheap)
# Visit: https://console.groq.com/keys

# 2. Set your API key
export GROQ_API_KEY="your-groq-api-key-here"

# 3. Install dependencies
poetry install

# 4. Run the summarizer
poetry run math-summarizer inputs/ch1.md --verbose
```

## 📋 What Changed

### New Dependencies Added
- `llama-index-llms-groq` - Groq API integration
- `llama-index-llms-together` - Together AI integration  
- `llama-index-llms-anthropic` - Anthropic Claude integration

### Configuration Changes
- **New provider system**: Choose between `groq`, `together`, `anthropic`, `huggingface`
- **Provider-specific API keys**: Each provider has its own key configuration
- **Default provider**: Now defaults to Groq for best cost/performance

### CLI Changes
- **New `--provider` option**: Explicitly choose your LLM provider
- **Updated `--model` option**: Now provider-specific (e.g., `llama-3.1-8b-instant` for Groq)
- **Better status reporting**: Shows which provider and model you're using

### Code Architecture Changes
- **Multi-provider summarizer**: `Summarizer` class now supports all providers
- **LlamaIndex integration**: Uses native LlamaIndex LLM classes for better performance
- **Fallback support**: Keeps original Hugging Face API as fallback option

## 🔧 Migration Steps

### If you were using Hugging Face API:

**Before:**
```bash
export HF_TOKEN="your-hf-token"
math-summarizer textbook.md --model facebook/bart-large-cnn
```

**After (Recommended - Groq):**
```bash
export GROQ_API_KEY="your-groq-key"
math-summarizer textbook.md --provider groq
```

**Or keep Hugging Face:**
```bash
export HF_TOKEN="your-hf-token"
math-summarizer textbook.md --provider huggingface --model facebook/bart-large-cnn
```

### Environment Variables

**Old `.env` format:**
```bash
HF_TOKEN=your-token
HF_MODEL_NAME=facebook/bart-large-cnn
```

**New `.env` format:**
```bash
LLM_PROVIDER=groq
GROQ_API_KEY=your-groq-key
GROQ_MODEL=llama-3.1-8b-instant
```

## 💰 Cost Benefits

| Provider | Cost vs OpenAI | Speed | Best For |
|----------|----------------|-------|----------|
| **Groq** | 90% cheaper | 10x faster | General use (recommended) |
| Together AI | 70% cheaper | 3x faster | Large models, technical content |
| Anthropic | Similar cost | Similar speed | Highest quality, complex reasoning |
| Hugging Face | Free tier | Variable | Experimentation |

## 🧪 Testing Your Setup

### 1. Test Configuration Loading
```bash
poetry run python -c "
from src.math_summarizer.config import Config
c = Config()
print(f'Provider: {c.llm_provider}')
print(f'Model: {c.groq_model if c.llm_provider == \"groq\" else \"N/A\"}')"
```

### 2. Test Different Providers
```bash
# Test Groq
export GROQ_API_KEY="your-key"
poetry run math-summarizer examples/sample_input.md --provider groq --verbose

# Test Together AI
export TOGETHER_API_KEY="your-key"
poetry run math-summarizer examples/sample_input.md --provider together --verbose

# Test Anthropic
export ANTHROPIC_API_KEY="your-key"
poetry run math-summarizer examples/sample_input.md --provider anthropic --verbose
```

### 3. Test Custom Models
```bash
# Use Groq's larger model
poetry run math-summarizer textbook.md --provider groq --model llama-3.1-70b-versatile

# Use Together AI's Mixtral
poetry run math-summarizer textbook.md --provider together --model mistralai/Mixtral-8x7B-Instruct-v0.1
```

## 🐛 Troubleshooting

### Import Errors
```bash
# Make sure all dependencies are installed
poetry install
```

### API Key Errors
```bash
# Check your environment variables
env | grep -E "(GROQ|TOGETHER|ANTHROPIC|HF)_"

# Or check your .env file
cat .env
```

### Model Not Found
```bash
# Use provider-specific models
# Groq: llama-3.1-8b-instant, llama-3.1-70b-versatile, mixtral-8x7b-32768
# Together: meta-llama/Llama-3-8b-chat-hf, mistralai/Mixtral-8x7B-Instruct-v0.1
# Anthropic: claude-3-haiku-20240307, claude-3-sonnet-20240229
```

### Performance Issues
```bash
# Use Groq for fastest processing
poetry run math-summarizer textbook.md --provider groq --model llama-3.1-8b-instant

# Reduce chunk size for large files
poetry run math-summarizer textbook.md --chunk-size 2000
```

## 📚 Documentation Updated

- ✅ `README.md` - Updated with multi-provider info
- ✅ `docs/usage.md` - Complete usage guide for all providers
- ✅ `API_SETUP.md` - Quick setup guide for each provider
- ✅ `.env` - Template with all provider configurations

## 🎯 Recommendations

1. **Start with Groq** - Best balance of speed, cost, and quality
2. **Use Together AI** for larger models when you need more capability
3. **Try Anthropic Claude** for highest quality on complex mathematical content
4. **Keep Hugging Face** as a fallback if other APIs have issues

The migration maintains backward compatibility while providing much better performance and cost options!
