"""Multi-provider LLM summarizer with clean architecture."""

from typing import List
import re

from ..config import Config
from ..chunking.chunker import Chunk
from ..utils.logger import setup_logger
from ..utils import LaTeXFixer, LaTeXDetector
from ..execution import ExecutionFactory

logger = setup_logger(__name__)

# Your system instructions
SYSTEM_INSTRUCTIONS = """
You are a mathematical content summarizer. Provide ONLY the final markdown summary - no reasoning, no thinking process, no explanations.

CRITICAL: Do not include any <think>, </think>, reasoning traces, or explanatory text. Output ONLY clean markdown.

Summarize the provided textbook chapter content using these rules:
1. Capture all theorems, propositions, definitions, lemmas, properties, conjectures in the provided chunk.
2. Do not include: exercises, proofs, notes, problems and sources, or any other text that is not a theorem, proposition, definition, lemma, property, or conjecture.
3. Use headings: "## [Section Number]. [Section Title]" for sections and "### [Section Number].[Subsection Number] [Subsection Title]" for subsections. 
4. Number all definitions and theorems using "###### Definition [Section].[Definition Number]" or "###### Theorem [Section].[Theorem Number]". 
5. For subsection items, include the subsection number, like "###### Definition 1.1.1". 
6. Use "$...$" for inline math and "$$...$$" for block equations. 
7. CRITICAL LaTeX RULES:
   - Each $$...$$ block must be on its own line with blank lines before and after
   - Never concatenate multiple $$...$$ blocks: $$eq1$$$$eq2$$ is WRONG
   - Always add spaces around inline math: $x$is wrong, $x$ is correct
   - Ensure every $$ has a matching closing $$
8. When referencing other definitions or theorems, use Obsidian's link format: "[[Title#Heading|Display Text]]".
9. Output clean markdown only - no meta-commentary, no reasoning process.

EXAMPLE OF EXPECTED OUTPUT:
## 3.1 Vector Spaces

###### Definition 3.1.1
Let $V$ be a non-empty set and $F$ be a field. Then $V$ is called a vector space over $F$ if it satisfies the following axioms:
1. **Closure under addition**: For all $\mathbf{u}, \mathbf{v} \in V$, we have $\mathbf{u} + \mathbf{v} \in V$
2. **Associativity**: $(\mathbf{u} + \mathbf{v}) + \mathbf{w} = \mathbf{u} + (\mathbf{v} + \mathbf{w})$ for all $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$
3. **Additive identity**: There exists $\mathbf{0} \in V$ such that $\mathbf{v} + \mathbf{0} = \mathbf{v}$ for all $\mathbf{v} \in V$
4. **Additive inverse**: For each $\mathbf{v} \in V$, there exists $-\mathbf{v} \in V$ such that $\mathbf{v} + (-\mathbf{v}) = \mathbf{0}$

###### Theorem 3.1.2
Let $V$ be a vector space over field $F$. Then:
$$\mathbf{0} \cdot \mathbf{v} = \mathbf{0}$$
for all $\mathbf{v} \in V$, where $\mathbf{0}$ on the left is the zero scalar and $\mathbf{0}$ on the right is the zero vector.

IMPORTANT: Each display math equation ($$...$$) should be on its own line with proper spacing. Never concatenate multiple $$...$$...$$...$$.

### 3.2 Linear Independence

###### Definition 3.2.1
A set of vectors $\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n\}$ in vector space $V$ is called **linearly independent** if the only solution to: $$c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_n\mathbf{v}_n = \mathbf{0}$$ is $c_1 = c_2 = \cdots = c_n = 0$.

---

Begin your response immediately with the markdown content following this exact format.
"""


class Summarizer:
    """Multi-provider LLM summarizer with clean architecture."""
    
    def __init__(self, config: Config):
        self.config = config
        self.latex_detector = LaTeXDetector()
        self.latex_fixer = LaTeXFixer(config)
        self.executor = ExecutionFactory.create_executor(config)
    
    
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
            # Use the executor for completion
            summary = self.executor.execute_completion(prompt)
            
            # Validate output and perform sanity checks
            summary = self._validate_output(summary)
            
            logger.debug(f"Generated summary for {chunk.chunk_id}: {len(summary)} characters")
            return summary
            
        except Exception as e:
            logger.error(f"Error summarizing chunk {chunk.chunk_id}: {str(e)}")
            # Return original content as fallback
            return f"## Error Processing Chunk {chunk.chunk_index}\n\n{chunk.content[:500]}..."
    
    
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
        
        # Sanity check 2: Detect and fix LaTeX formatting issues
        latex_issues = self.latex_detector.detect_all_issues(text)
        
        if latex_issues:
            issue_summary = self.latex_detector.get_issue_summary(latex_issues)
            logger.info(f"Detected {issue_summary['total_issues']} LaTeX issues: {issue_summary['issue_types']}")
            
            # Fix the issues using LLM
            text = self.latex_fixer.fix_latex(text, latex_issues)
        else:
            logger.debug("No LaTeX issues detected")
        
        # Sanity check 3: Check for mathematical content indicators
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
    