"""LaTeX formatting repair utility using LLM."""

import re
import logging
from typing import Optional, List
from ..llm_clients import LLMClientFactory, BaseLLMClient
from ..config import Config
from .latex_detector import LaTeXIssue

logger = logging.getLogger(__name__)


class LaTeXFixer:
    """A utility class for fixing LaTeX formatting issues using an LLM."""
    
    def __init__(self, config: Optional[Config] = None):
        """Initialize the LaTeX fixer with LLM client."""
        if config is None:
            from ..config import config as default_config
            config = default_config
        
        self.config = config
        self._client = None
        
        # LaTeX fixing prompt template
        self.fix_prompt = """You are a LaTeX formatting specialist. Fix any LaTeX formatting issues in the following mathematical content.

CRITICAL RULES:
1. Output ONLY the corrected content - no explanations, no reasoning, no additional text
2. Fix common LaTeX issues:
   - Unmatched $$ delimiters (ensure every $$ has a closing $$)
   - Concatenated display math blocks: $$eq1$$eq2$$ → $$eq1$$\n\n$$eq2$$
   - Missing spaces around inline math: $x$is → $x$ is
   - Empty math environments: $$ $$ → (remove)
   - Excessive newlines: \n\n\n\n → \n\n
3. Remove any reasoning traces like "Let me...", "I need to...", <think>...</think>
4. Preserve all mathematical content and structure exactly
5. Each $$...$$ block should be on its own line with proper spacing

CONTENT TO FIX:
{content}

OUTPUT THE CORRECTED CONTENT:"""
    
    def _get_client(self) -> BaseLLMClient:
        """Get or create the LLM client for LaTeX fixing."""
        if self._client is None:
            # Use the main client for LaTeX fixing
            try:
                self._client = LLMClientFactory.create_main_client(self.config)
                logger.info(f"LaTeX fixer initialized with: {self._client.get_config_summary()}")
            except Exception as e:
                logger.error(f"Failed to initialize LaTeX fixer client: {e}")
                raise
        
        return self._client
    
    def fix_latex(self, text: str, issues: Optional[List[LaTeXIssue]] = None) -> str:
        """
        Fix LaTeX formatting issues in the text using an LLM.
        
        Args:
            text: The text content to fix
            issues: Optional list of detected issues to guide fixing
            
        Returns:
            The fixed text content
        """
        if issues is None or len(issues) == 0:
            logger.debug("No LaTeX issues provided, returning original text")
            return text
        
        # Create a detailed issue description for the LLM
        issue_descriptions = []
        for issue in issues:
            desc = f"- {issue.issue_type}: {issue.description}"
            if issue.line_number:
                desc += f" (line {issue.line_number})"
            if issue.context:
                desc += f" - Context: {issue.context}"
            issue_descriptions.append(desc)
        
        issues_text = "\n".join(issue_descriptions)
        
        logger.info(f"Fixing {len(issues)} LaTeX issues with LLM")
        
        try:
            client = self._get_client()
            
            # Enhanced prompt with specific issue information
            enhanced_prompt = f"""You are a LaTeX formatting specialist. Fix the specific LaTeX formatting issues listed below in the mathematical content.

DETECTED ISSUES TO FIX:
{issues_text}

CRITICAL RULES:
1. Output ONLY the corrected content - no explanations, no reasoning, no additional text
2. Fix ONLY the issues listed above
3. Preserve all mathematical content and structure exactly

CONTENT TO FIX:
{text}

OUTPUT THE CORRECTED CONTENT:"""
            
            fixed_text = client.complete(enhanced_prompt).strip()
            
            # Basic validation
            if len(fixed_text) < len(text) * 0.3:
                logger.warning("LaTeX fixer returned suspiciously short text, using original")
                return text
            
            logger.info("LaTeX formatting completed successfully")
            return fixed_text
            
        except Exception as e:
            logger.error(f"LaTeX fixing failed: {e}")
            logger.info("Returning original text due to LaTeX fixer failure")
            return text