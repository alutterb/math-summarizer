"""Hugging Face LLM client implementation."""

from .base import BaseLLMClient
import logging

logger = logging.getLogger(__name__)


class HuggingFaceClient(BaseLLMClient):
    """Hugging Face Inference API client implementation."""
    
    def _initialize_client(self) -> None:
        """Initialize the Hugging Face client."""
        try:
            from huggingface_hub import InferenceClient
            
            self._client = InferenceClient(
                model=self.model,
                token=self.api_key,
                timeout=self.timeout
            )
            logger.info(f"Hugging Face client initialized with model: {self.model}")
            
        except ImportError:
            raise ImportError("huggingface_hub is required for Hugging Face support. Install with: pip install huggingface_hub")
        except Exception as e:
            logger.error(f"Failed to initialize Hugging Face client: {str(e)}")
            raise
    
    def _generate_text(self, prompt: str) -> str:
        """Generate text using Hugging Face Inference API."""
        response = self._client.text_generation(
            prompt,
            max_new_tokens=self.max_tokens,
            temperature=self.temperature,
            do_sample=self.temperature > 0,
            return_full_text=False,
        )
        return response
