"""Local execution for running models locally."""

import logging
from .base import BaseExecutor
from ..config import Config

logger = logging.getLogger(__name__)


class LocalExecutor(BaseExecutor):
    """Handles execution using local model instances."""
    
    def __init__(self, config: Config):
        super().__init__(config)
        self.llm = None
    
    def initialize(self) -> None:
        """Initialize local model."""
        logger.info(f"Initializing local mode with model: {self.config.hf_model_name}")
        
        try:
            # Check if local dependencies are available
            import torch
            from llama_index.llms.huggingface import HuggingFaceLLM
            from llama_index.core import Settings
            
            # Determine device
            if self.config.device == "auto":
                device = "cuda" if torch.cuda.is_available() else "cpu"
            else:
                device = self.config.device
            
            logger.info(f"Using device: {device}")
            
            # Initialize local Hugging Face LLM
            self.llm = HuggingFaceLLM(
                model_name=self.config.hf_model_name,
                tokenizer_name=self.config.hf_model_name,
                context_window=4096,
                max_new_tokens=self.config.max_new_tokens,
                model_kwargs={
                    "torch_dtype": torch.float16,
                    "load_in_8bit": self.config.load_in_8bit,
                    "load_in_4bit": self.config.load_in_4bit,
                    "device_map": "auto" if device == "cuda" else None,
                },
                generate_kwargs={
                    "temperature": self.config.temperature,
                    "do_sample": self.config.do_sample,
                    "top_p": self.config.top_p,
                    "top_k": self.config.top_k,
                    "pad_token_id": 50256,
                },
                device_map="auto" if device == "cuda" else None,
                tokenizer_kwargs={"max_length": 4096},
            )
            
            Settings.llm = self.llm
            logger.info("Successfully initialized local Hugging Face model")
            
        except ImportError as e:
            logger.error("Local mode requires additional dependencies. Install with: poetry install -E local")
            raise ImportError("Missing local execution dependencies. Run: poetry install -E local") from e
        except Exception as e:
            logger.error(f"Failed to initialize local model: {str(e)}")
            raise
    
    def execute_completion(self, prompt: str) -> str:
        """Execute completion using the local model."""
        if self.llm is None:
            raise RuntimeError("Local executor not initialized. Call initialize() first.")
        
        response = self.llm.complete(prompt)
        return response.text.strip()
    
    def get_execution_info(self) -> dict:
        """Get local execution information."""
        if self.llm is None:
            return {"mode": "local", "status": "not_initialized"}
        
        return {
            "mode": "local",
            "model": self.config.hf_model_name,
            "device": self.config.device,
            "status": "initialized"
        }
