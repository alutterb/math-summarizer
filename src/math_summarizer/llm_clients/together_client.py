"""Together AI LLM client implementation."""

from .base import BaseLLMClient
import logging

logger = logging.getLogger(__name__)


class TogetherClient(BaseLLMClient):
    """Together AI client implementation."""
    
    def _initialize_client(self) -> None:
        """Initialize the Together AI client."""
        try:
            from llama_index.llms.together import TogetherLLM
            
            self._client = TogetherLLM(
                model=self.model,
                api_key=self.api_key,
                temperature=self.temperature,
                max_tokens=self.max_tokens,
            )
            logger.info(f"Together AI client initialized with model: {self.model}")
            
        except ImportError:
            raise ImportError("llama-index-llms-together is required for Together AI support. Install with: pip install llama-index-llms-together")
        except Exception as e:
            logger.error(f"Failed to initialize Together AI client: {str(e)}")
            raise
    
    def _generate_text(self, prompt: str) -> str:
        """Generate text using Together AI API."""
        response = self._client.complete(prompt)
        return response.text
