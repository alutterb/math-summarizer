"""Configuration management for the Math Summarizer."""

import os
from typing import Optional
from pydantic_settings import BaseSettings
from pydantic import Field
from dotenv import load_dotenv

load_dotenv()


class Config(BaseSettings):
    """Configuration settings for the Math Summarizer."""
    
    # LLM Provider Configuration
    llm_provider: str = Field(default="groq", env="LLM_PROVIDER")  # groq, together, anthropic, huggingface
    
    # Hugging Face Model Configuration (when provider=huggingface)
    hf_model_name: str = Field(default="microsoft/DialoGPT-medium", env="HF_MODEL_NAME")
    hf_token: Optional[str] = Field(default=None, env="HF_TOKEN")
    
    # Groq Configuration (when provider=groq)
    groq_api_key: Optional[str] = Field(default=None, env="GROQ_API_KEY")
    groq_model: str = Field(default="llama-3.1-8b-instant", env="GROQ_MODEL")
    
    # Together AI Configuration (when provider=together)
    together_api_key: Optional[str] = Field(default=None, env="TOGETHER_API_KEY")
    together_model: str = Field(default="meta-llama/Llama-3-8b-chat-hf", env="TOGETHER_MODEL")
    
    # Anthropic Configuration (when provider=anthropic)
    anthropic_api_key: Optional[str] = Field(default=None, env="ANTHROPIC_API_KEY")
    anthropic_model: str = Field(default="claude-3-haiku-20240307", env="ANTHROPIC_MODEL")
    
    # LaTeX Cleanup Agent Configuration (separate lightweight model)
    latex_cleanup_provider: str = Field(default="groq", env="LATEX_CLEANUP_PROVIDER")
    latex_cleanup_model: str = Field(default="llama-3.1-8b-instant", env="LATEX_CLEANUP_MODEL")
    latex_cleanup_api_key: Optional[str] = Field(default=None, env="LATEX_CLEANUP_API_KEY")  # Falls back to main provider key if not set
    
    # Execution Mode: "api" for remote APIs, "local" for local execution
    execution_mode: str = Field(default="api", env="EXECUTION_MODE")
    
    # API Configuration (for remote execution)
    api_timeout: int = Field(default=120, env="API_TIMEOUT")  # seconds
    max_retries: int = Field(default=3, env="MAX_RETRIES")
    
    # Local Configuration (only used when execution_mode="local")
    device: str = Field(default="auto", env="DEVICE")  # auto, cpu, cuda
    load_in_8bit: bool = Field(default=False, env="LOAD_IN_8BIT")
    load_in_4bit: bool = Field(default=False, env="LOAD_IN_4BIT")
    
    # Generation Configuration (works for both modes)
    max_new_tokens: int = Field(default=512, env="MAX_NEW_TOKENS")
    
    # Chunking Configuration
    max_chunk_size: int = Field(default=4000, env="MAX_CHUNK_SIZE")
    chunk_overlap: int = Field(default=200, env="CHUNK_OVERLAP")
    min_chunk_size: int = Field(default=100, env="MIN_CHUNK_SIZE")
    
    # Generation Configuration  
    temperature: float = Field(default=0.7, env="TEMPERATURE")
    do_sample: bool = Field(default=True, env="DO_SAMPLE")
    top_p: float = Field(default=0.9, env="TOP_P")
    top_k: int = Field(default=50, env="TOP_K")
    
    # Output Configuration
    output_format: str = Field(default="markdown", env="OUTPUT_FORMAT")
    preserve_math_notation: bool = Field(default=True, env="PRESERVE_MATH_NOTATION")
    
    # Logging
    log_level: str = Field(default="INFO", env="LOG_LEVEL")
    log_file: Optional[str] = Field(default="logs/math_summarizer.log", env="LOG_FILE")
    
    class Config:
        env_file = ".env"
        case_sensitive = False


# Global configuration instance
config = Config()
