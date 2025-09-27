"""Factory for creating LLM clients."""

from typing import Optional, Dict, Any
import logging

from .base import BaseLLMClient
from .groq_client import GroqClient
from .together_client import TogetherClient
from .anthropic_client import AnthropicClient
from .huggingface_client import HuggingFaceClient
from ..config import Config

logger = logging.getLogger(__name__)


class LLMClientFactory:
    """Factory for creating LLM clients based on provider configuration."""
    
    _CLIENT_CLASSES = {
        'groq': GroqClient,
        'together': TogetherClient,
        'anthropic': AnthropicClient,
        'huggingface': HuggingFaceClient,
    }
    
    @classmethod
    def create_main_client(cls, config: Config) -> BaseLLMClient:
        """Create the main summarization client."""
        provider = config.llm_provider.lower()
        
        # Get provider-specific configuration
        model, api_key = cls._get_main_provider_config(config, provider)
        
        return cls._create_client(
            provider=provider,
            model=model,
            api_key=api_key,
            temperature=config.temperature,
            max_tokens=config.max_new_tokens,
            timeout=config.api_timeout,
            max_retries=config.max_retries
        )
    
    
    @classmethod
    def _create_client(cls, provider: str, model: str, api_key: str, **kwargs) -> BaseLLMClient:
        """Create a client instance for the specified provider."""
        if provider not in cls._CLIENT_CLASSES:
            raise ValueError(f"Unsupported provider: {provider}. Supported providers: {list(cls._CLIENT_CLASSES.keys())}")
        
        client_class = cls._CLIENT_CLASSES[provider]
        client = client_class(model=model, api_key=api_key, **kwargs)
        
        logger.info(f"Created {provider} client: {client.get_config_summary()}")
        return client
    
    @classmethod
    def _get_main_provider_config(cls, config: Config, provider: str) -> tuple[str, str]:
        """Get model and API key for main provider."""
        provider_configs = {
            'groq': (config.groq_model, config.groq_api_key),
            'together': (config.together_model, config.together_api_key),
            'anthropic': (config.anthropic_model, config.anthropic_api_key),
            'huggingface': (config.hf_model_name, config.hf_token),
        }
        
        if provider not in provider_configs:
            raise ValueError(f"Unknown provider: {provider}")
        
        model, api_key = provider_configs[provider]
        if not api_key:
            raise ValueError(f"API key not configured for {provider} provider")
        
        return model, api_key
    
    
    @classmethod
    def get_supported_providers(cls) -> list[str]:
        """Get list of supported providers."""
        return list(cls._CLIENT_CLASSES.keys())
