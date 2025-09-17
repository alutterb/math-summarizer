"""Multi-provider LLM summarizer with clean architecture."""

from typing import List, Optional
import re
import time
from llama_index.core import Settings

from ..config import Config
from ..chunking.chunker import Chunk
from ..utils.logger import setup_logger
from ..llm_clients import LLMClientFactory, BaseLLMClient

logger = setup_logger(__name__)

# Your system instructions
SYSTEM_INSTRUCTIONS = """
You are a mathematical content summarizer. Provide ONLY the final markdown summary - no reasoning, no thinking process, no explanations.

CRITICAL: Do not include any <think>, </think>, reasoning traces, or explanatory text. Output ONLY clean markdown.

Summarize the provided textbook chapter content using these rules:
1. Capture all theorems, propositions, definitions, lemmas, properties, conjectures in the provided chunk.
2. Use headings: "## [Section Number]. [Section Title]" for sections and "### [Section Number].[Subsection Number] [Subsection Title]" for subsections. 
3. Within each section, include "##### Summary" followed by a brief summary, and "##### Keywords" followed by relevant keywords separated by commas. 
4. Number all definitions and theorems using "###### Definition [Section].[Definition Number]" or "###### Theorem [Section].[Theorem Number]". 
5. For subsection items, include the subsection number, like "###### Definition 1.1.1". 
6. Use "$...$" for inline math and "$$...$$" for block equations. 
7. When referencing definitions or theorems, use Obsidian's link format: "[[Title#Heading|Display Text]]".
8. Output clean markdown only - no meta-commentary, no reasoning process.

Begin your response immediately with the markdown content.
"""


class Summarizer:
    """Multi-provider LLM summarizer with clean architecture."""
    
    def __init__(self, config: Config):
        self.config = config
        
        if config.execution_mode.lower() == "api":
            self._init_api_mode()
        else:
            self._init_local_mode()
    
    def _init_api_mode(self):
        """Initialize API mode using factory pattern."""
        logger.info(f"Initializing API mode with provider: {self.config.llm_provider}")
        
        try:
            # Create main client
            self.main_client = LLMClientFactory.create_main_client(self.config)
            logger.info(f"Main client initialized: {self.main_client.get_config_summary()}")
            
            # Create LaTeX cleanup client
            self.latex_client = LLMClientFactory.create_latex_cleanup_client(self.config)
            if self.latex_client:
                logger.info(f"LaTeX cleanup client initialized: {self.latex_client.get_config_summary()}")
            else:
                logger.info("No dedicated LaTeX cleanup client - will use main client as fallback")
            
        except Exception as e:
            logger.error(f"Failed to initialize API clients: {str(e)}")
            raise
    
    
    def _init_local_mode(self):
        """Initialize for local execution."""
        logger.info(f"Initializing local mode with model: {self.config.hf_model_name}")
        
        try:
            # Check if local dependencies are available
            import torch
            from llama_index.llms.huggingface import HuggingFaceLLM
            
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
    
    def summarize_chunk(self, chunk: Chunk) -> str:
        """Summarize a single chunk using the configured LLM provider."""
        logger.debug(f"Summarizing chunk: {chunk.chunk_id} ({chunk.token_count} tokens)")
        
        prompt = f"""
{SYSTEM_INSTRUCTIONS}

Content to process:
{chunk.content}

Please process this content following the rules above:
"""
        
        try:
            if self.config.execution_mode.lower() == "api":
                # Use the main client for summarization
                summary = self.main_client.complete(prompt)
            else:
                summary = self._summarize_with_local(prompt)
            
            # Validate output and perform sanity checks
            summary = self._validate_output(summary)
            
            logger.debug(f"Generated summary for {chunk.chunk_id}: {len(summary)} characters")
            return summary
            
        except Exception as e:
            logger.error(f"Error summarizing chunk {chunk.chunk_id}: {str(e)}")
            # Return original content as fallback
            return f"## Error Processing Chunk {chunk.chunk_index}\n\n{chunk.content[:500]}..."
    
    def _summarize_with_llamaindex(self, prompt: str) -> str:
        """Summarize using LlamaIndex LLM (Groq, Together, Anthropic)."""
        response = self.llm.complete(prompt)
        return response.text.strip()
    
    def _summarize_with_hf_api(self, prompt: str) -> str:
        """Summarize using Hugging Face Inference API (legacy method)."""
        for attempt in range(self.config.max_retries):
            try:
                logger.debug(f"API attempt {attempt + 1}/{self.config.max_retries}")
                
                # Use text generation
                response = self.client.text_generation(
                    prompt,
                    max_new_tokens=self.config.max_new_tokens,
                    temperature=self.config.temperature,
                    do_sample=self.config.do_sample,
                    top_p=self.config.top_p,
                    top_k=self.config.top_k,
                    return_full_text=False,
                )
                
                return response
                
            except Exception as e:
                logger.warning(f"API attempt {attempt + 1} failed: {str(e)}")
                if attempt < self.config.max_retries - 1:
                    wait_time = 2 ** attempt  # Exponential backoff
                    logger.info(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    raise
    
    def _summarize_with_local(self, prompt: str) -> str:
        """Summarize using local model."""
        response = self.llm.complete(prompt)
        return response.text.strip()
    
    def summarize_chunks(self, chunks: List[Chunk]) -> List[str]:
        """Summarize multiple chunks."""
        logger.info(f"Summarizing {len(chunks)} chunks using {self.config.execution_mode} mode")
        
        summaries = []
        for i, chunk in enumerate(chunks, 1):
            logger.info(f"Processing chunk {i}/{len(chunks)}")
            summary = self.summarize_chunk(chunk)
            summaries.append(summary)
        
        logger.info(f"Completed summarization of {len(chunks)} chunks")
        return summaries
    
    def _validate_output(self, text: str) -> str:
        """Validate summarizer output and perform sanity checks."""
        import re
        
        # Sanity check 1: Check if response has content
        if not text or len(text.strip()) == 0:
            logger.warning("Missing content warning: Summarizer returned empty response")
            return "## Error: Empty Response\n\nThe summarizer returned no content for this chunk."
        
        # Sanity check 2: Remove reasoning traces and unwanted content
        original_length = len(text)
        
        # Remove <think> blocks completely
        text = re.sub(r'<think>.*?</think>', '', text, flags=re.DOTALL)
        
        # Remove other reasoning patterns
        reasoning_patterns = [
            r'Let me.*?(?=\n#|\n\n|$)',  # "Let me think about this..."
            r'I need to.*?(?=\n#|\n\n|$)',  # "I need to analyze..."
            r'First.*?(?=\n#|\n\n|$)',  # "First, let me..."
            r'Okay.*?(?=\n#|\n\n|$)',   # "Okay, so I need to..."
            r'Now.*?(?=\n#|\n\n|$)',    # "Now I'll process..."
        ]
        
        for pattern in reasoning_patterns:
            text = re.sub(pattern, '', text, flags=re.DOTALL)
        
        # Clean up multiple newlines
        text = re.sub(r'\n{3,}', '\n\n', text)
        text = text.strip()
        
        if len(text) < original_length * 0.5:
            logger.warning(f"Cleaned reasoning traces: reduced from {original_length} to {len(text)} characters")
        
        # Sanity check 3: Check for LaTeX formatting issues
        latex_issues = self._detect_latex_issues(text)
        
        if latex_issues:
            logger.warning(f"LaTeX formatting issues detected: {latex_issues}")
            text = self._fix_latex_with_agent(text, latex_issues)
        
        # Sanity check 4: Check for mathematical content indicators
        mathematical_terms = ['proposition', 'theorem', 'definition', 'corollary', 'lemma']
        text_lower = text.lower()
        
        found_terms = [term for term in mathematical_terms if term in text_lower]
        
        if not found_terms:
            logger.warning(f"Format warning: Summary may be missing mathematical content. "
                         f"Expected to find at least one of: {', '.join(mathematical_terms)}")
            # Add a warning to the output but don't fail
            warning_note = f"\n\n> **⚠️ Warning**: This summary may be incomplete. Expected mathematical content (definitions, theorems, propositions, corollaries) not detected.\n\n"
            text = warning_note + text
        else:
            logger.debug(f"Sanity check passed: Found mathematical terms: {', '.join(found_terms)}")
        
        return text
    
    def _detect_latex_issues(self, text: str) -> list:
        """Detect common LaTeX formatting issues."""
        import re
        issues = []
        
        # Issue 1: Spaces inside inline math delimiters ($ x = 5 $ instead of $x = 5$)
        if re.search(r'\$\s+[^$]+\s+\$', text):
            issues.append("spaces_in_inline_math")
        
        # Issue 2: Spaces inside display math delimiters ($$ x = 5 $$ instead of $$x = 5$$)
        if re.search(r'\$\$\s+[^$]+\s+\$\$', text):
            issues.append("spaces_in_display_math")
        
        # Issue 3: Unmatched dollar signs
        dollar_count = text.count('$')
        if dollar_count % 2 != 0:
            issues.append("unmatched_dollar_signs")
        
        # Issue 4: Empty math environments
        if re.search(r'\$\s*\$', text) or re.search(r'\$\$\s*\$\$', text):
            issues.append("empty_math_environments")
        
        # Issue 5: Broken LaTeX commands (missing backslashes, etc.)
        common_commands = ['frac', 'sqrt', 'sum', 'int', 'lim', 'nabla', 'alpha', 'beta', 'gamma', 'theta', 'lambda', 'mu', 'sigma', 'pi']
        for cmd in common_commands:
            if re.search(rf'\b{cmd}\b(?![a-zA-Z])', text) and not re.search(rf'\\{cmd}', text):
                issues.append(f"missing_backslash_{cmd}")
        
        # Issue 6: Incorrect matrix/equation formatting
        if re.search(r'\$[^$]*\n[^$]*\$', text, re.MULTILINE):
            issues.append("multiline_inline_math")
        
        return issues
    
    def _fix_latex_with_agent(self, text: str, issues: list) -> str:
        """Use a specialized agent to fix LaTeX formatting issues."""
        logger.info(f"Calling LaTeX cleanup agent for issues: {', '.join(issues)}")
        
        cleanup_prompt = f"""You are a LaTeX formatting specialist. Fix the LaTeX formatting issues in the following mathematical content.

DETECTED ISSUES: {', '.join(issues)}

RULES FOR FIXING:
MAIN RULE: Do not add explanations,reasoning or ANY additional text beyond what was corrected - output only the corrected content.
1. Fix spaces in math delimiters: "$ x = 5 $" → "$x = 5$"
2. Fix spaces in display math: "$$ equation $$" → "$$equation$$"
3. Match unmatched dollar signs
4. Remove empty math environments: "$$ $$" → ""
5. Add missing backslashes to LaTeX commands: "alpha" → "\\alpha"
6. Convert multiline inline math to display math when appropriate
7. Preserve all non-mathematical content exactly as is

CONTENT TO FIX:
{text}

OUTPUT THE CORRECTED CONTENT:"""

        try:
            if self.latex_client is not None:
                # Use dedicated LaTeX cleanup client
                logger.info(f"Using dedicated LaTeX cleanup client: {self.latex_client.get_config_summary()}")
                cleaned_text = self.latex_client.complete(cleanup_prompt)
            else:
                # Fallback to main client
                logger.info("Using main client for LaTeX cleanup (no dedicated client)")
                cleaned_text = self.main_client.complete(cleanup_prompt)
            
            logger.info("LaTeX cleanup completed successfully")
            return cleaned_text.strip()
            
        except Exception as e:
            logger.error(f"LaTeX cleanup failed: {str(e)}")
            logger.info("Applying basic LaTeX fixes as fallback")
            return self._apply_basic_latex_fixes(text)
    
    
    def _apply_basic_latex_fixes(self, text: str) -> str:
        """Apply basic LaTeX fixes as fallback when agent fails."""
        import re
        
        # Fix spaces in inline math
        text = re.sub(r'\$\s+([^$]+?)\s+\$', r'$\1$', text)
        
        # Fix spaces in display math
        text = re.sub(r'\$\$\s+([^$]+?)\s+\$\$', r'$$\1$$', text)
        
        # Remove empty math environments
        text = re.sub(r'\$\s*\$', '', text)
        text = re.sub(r'\$\$\s*\$\$', '', text)
        
        # Add backslashes to common math commands (basic patterns)
        common_fixes = {
            r'\balpha\b': r'\\alpha',
            r'\bbeta\b': r'\\beta',
            r'\bgamma\b': r'\\gamma',
            r'\btheta\b': r'\\theta',
            r'\blambda\b': r'\\lambda',
            r'\bmu\b': r'\\mu',
            r'\bsigma\b': r'\\sigma',
            r'\bpi\b': r'\\pi',
            r'\bnabla\b': r'\\nabla'
        }
        
        for pattern, replacement in common_fixes.items():
            # Only fix if it's inside math delimiters
            def fix_in_math(match):
                content = match.group(0)
                for p, r in common_fixes.items():
                    content = re.sub(p, r, content)
                return content
            
            text = re.sub(r'\$[^$]*\$', fix_in_math, text)
            text = re.sub(r'\$\$[^$]*\$\$', fix_in_math, text)
        
        logger.info("Applied basic LaTeX fixes as fallback")
        return text