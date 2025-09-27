"""LaTeX issue detection utility."""

import re
import logging
from typing import List, Dict, Any
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class LaTeXIssue:
    """Represents a detected LaTeX formatting issue."""
    issue_type: str
    description: str
    line_number: int = None
    context: str = None
    severity: str = "medium"  # low, medium, high


class LaTeXDetector:
    """A utility class for detecting LaTeX formatting issues."""
    
    def __init__(self):
        """Initialize the LaTeX detector with issue patterns."""
        self.common_commands = [
            'frac', 'sqrt', 'sum', 'int', 'lim', 'nabla', 
            'alpha', 'beta', 'gamma', 'delta', 'theta', 'lambda', 
            'mu', 'sigma', 'pi', 'phi', 'psi', 'omega', 'infty',
            'partial', 'times', 'cdot', 'ldots', 'geq', 'leq'
        ]
        
        self.common_environments = [
            'equation', 'align', 'gather', 'matrix', 'bmatrix', 'pmatrix',
            'cases', 'split', 'multline', 'eqnarray', 'aligned'
        ]
    
    def detect_all_issues(self, text: str) -> List[LaTeXIssue]:
        """
        Detect all LaTeX formatting issues in the text.
        
        Returns:
            List of LaTeXIssue objects describing each problem found
        """
        issues = []
        lines = text.split('\n')
        
        # Issue detection methods
        issues.extend(self._detect_math_delimiter_issues(text, lines))
        issues.extend(self._detect_environment_issues(text, lines))
        issues.extend(self._detect_command_issues(text, lines))
        issues.extend(self._detect_formatting_issues(text, lines))
        issues.extend(self._detect_content_issues(text, lines))
        
        logger.info(f"Detected {len(issues)} LaTeX issues")
        return issues
    
    def get_issue_summary(self, issues: List[LaTeXIssue]) -> Dict[str, Any]:
        """Get a summary of detected issues by type and severity."""
        summary = {
            "total_issues": len(issues),
            "by_type": {},
            "by_severity": {"low": 0, "medium": 0, "high": 0},
            "issue_types": []
        }
        
        for issue in issues:
            # Count by type
            if issue.issue_type not in summary["by_type"]:
                summary["by_type"][issue.issue_type] = 0
            summary["by_type"][issue.issue_type] += 1
            
            # Count by severity
            summary["by_severity"][issue.severity] += 1
            
            # Collect unique issue types
            if issue.issue_type not in summary["issue_types"]:
                summary["issue_types"].append(issue.issue_type)
        
        return summary
    
    def _detect_math_delimiter_issues(self, text: str, lines: List[str]) -> List[LaTeXIssue]:
        """Detect issues with math delimiters ($, $$)."""
        issues = []
        
        # Issue 1: Unmatched dollar signs
        dollar_count = text.count('$')
        if dollar_count % 2 != 0:
            issues.append(LaTeXIssue(
                issue_type="unmatched_dollar_signs",
                description=f"Unmatched dollar signs (found {dollar_count}, should be even)",
                severity="high"
            ))
        
        # Issue 2: Concatenated display math blocks (more specific pattern)
        concatenated_matches = re.finditer(r'\$\$[^$]+\$\$\$\$[^$]+\$\$', text)
        for match in concatenated_matches:
            line_num = text[:match.start()].count('\n') + 1
            issues.append(LaTeXIssue(
                issue_type="concatenated_display_math",
                description="Multiple display math blocks concatenated without spacing",
                line_number=line_num,
                context=match.group()[:100] + "..." if len(match.group()) > 100 else match.group(),
                severity="high"
            ))
        
        # Issue 3: Display math blocks not on their own line (only flag very long equations)
        # Short display math like $$V(x^*)$$ is often acceptable inline
        inline_display_matches = re.finditer(r'[^\n]\$\$[^$]{20,}\$\$[^\n]', text)
        for match in inline_display_matches:
            line_num = text[:match.start()].count('\n') + 1
            issues.append(LaTeXIssue(
                issue_type="inline_display_math",
                description="Long display math block should be on its own line",
                line_number=line_num,
                context=match.group()[:50] + "..." if len(match.group()) > 50 else match.group(),
                severity="medium"
            ))
        
        # Issue 4: Excessive spaces in math delimiters
        for i, line in enumerate(lines, 1):
            if re.search(r'\$\s{2,}[^$]+\s{2,}\$', line):
                issues.append(LaTeXIssue(
                    issue_type="excessive_spaces_inline_math",
                    description="Excessive spaces inside inline math delimiters",
                    line_number=i,
                    context=line.strip(),
                    severity="low"
                ))
            
            if re.search(r'\$\$\s{2,}[^$]+\s{2,}\$\$', line):
                issues.append(LaTeXIssue(
                    issue_type="excessive_spaces_display_math",
                    description="Excessive spaces inside display math delimiters",
                    line_number=i,
                    context=line.strip(),
                    severity="low"
                ))
        
        # Issue 5: Empty math environments (only whitespace-only)
        empty_math_matches = re.finditer(r'\$\s+\$|\$\$\s+\$\$', text)
        for match in empty_math_matches:
            line_num = text[:match.start()].count('\n') + 1
            issues.append(LaTeXIssue(
                issue_type="empty_math_environment",
                description="Empty math environment",
                line_number=line_num,
                context=match.group(),
                severity="medium"
            ))
        
        # Issue 6: Missing spaces around inline math (only obvious cases)
        # Only flag cases where there's clearly no space on BOTH sides
        missing_space_matches = re.finditer(r'[a-zA-Z]\$[^$]+\$[a-zA-Z]', text)
        for match in missing_space_matches:
            # Additional check: make sure it's not a valid pattern like "where$V$is"
            context = match.group()
            # Only flag if the math is a single variable without spaces on both sides
            if re.match(r'[a-zA-Z]\$[a-zA-Z_]\$[a-zA-Z]', context):
                line_num = text[:match.start()].count('\n') + 1
                issues.append(LaTeXIssue(
                    issue_type="missing_spaces_around_math",
                    description="Missing spaces around inline math",
                    line_number=line_num,
                    context=context,
                    severity="medium"
                ))
        
        return issues
    
    def _detect_environment_issues(self, text: str, lines: List[str]) -> List[LaTeXIssue]:
        """Detect issues with LaTeX environments."""
        issues = []
        
        # Check for unmatched \begin{} and \end{} blocks
        for env in self.common_environments:
            begin_pattern = rf'\\begin\{{{env}\}}'
            end_pattern = rf'\\end\{{{env}\}}'
            
            begin_matches = list(re.finditer(begin_pattern, text))
            end_matches = list(re.finditer(end_pattern, text))
            
            begin_count = len(begin_matches)
            end_count = len(end_matches)
            
            if begin_count != end_count:
                issues.append(LaTeXIssue(
                    issue_type="unmatched_environment",
                    description=f"Unmatched \\begin{{{env}}} and \\end{{{env}}} blocks ({begin_count} begins, {end_count} ends)",
                    severity="high"
                ))
        
        return issues
    
    def _detect_command_issues(self, text: str, lines: List[str]) -> List[LaTeXIssue]:
        """Detect issues with LaTeX commands."""
        issues = []
        
        # Check for missing backslashes on common commands (only inside math)
        math_blocks = re.findall(r'\$[^$]+\$|\$\$[^$]+\$\$', text)
        for block in math_blocks:
            for cmd in self.common_commands:
                if re.search(rf'\b{cmd}\b(?![a-zA-Z])', block) and not re.search(rf'\\{cmd}', block):
                    issues.append(LaTeXIssue(
                        issue_type="missing_backslash",
                        description=f"Missing backslash for command '{cmd}' in math block",
                        context=block[:50] + "..." if len(block) > 50 else block,
                        severity="medium"
                    ))
        
        return issues
    
    def _detect_formatting_issues(self, text: str, lines: List[str]) -> List[LaTeXIssue]:
        """Detect general formatting issues."""
        issues = []
        
        # Issue 1: Excessive newlines (4 or more consecutive)
        excessive_newlines = re.finditer(r'\n{4,}', text)
        for match in excessive_newlines:
            line_num = text[:match.start()].count('\n') + 1
            newline_count = len(match.group())
            issues.append(LaTeXIssue(
                issue_type="excessive_newlines",
                description=f"Excessive newlines ({newline_count} consecutive)",
                line_number=line_num,
                severity="low"
            ))
        
        # Issue 2: Very obvious missing spaces around operators (be conservative)
        for i, line in enumerate(lines, 1):
            # Only check for very obvious cases like "x=y" without any spaces
            if re.search(r'\$[^$]*[a-zA-Z0-9]=[a-zA-Z0-9][^$]*\$', line):
                # Additional check: make sure it's not something like "i=1" which is often OK
                if not re.search(r'[a-zA-Z]=\d+|\d+=[a-zA-Z]', line):
                    issues.append(LaTeXIssue(
                        issue_type="missing_operator_spacing",
                        description="Missing spaces around equals sign in math",
                        line_number=i,
                        context=line.strip(),
                        severity="low"
                    ))
        
        return issues
    
    def _detect_content_issues(self, text: str, lines: List[str]) -> List[LaTeXIssue]:
        """Detect content-related issues (reasoning traces, etc.)."""
        issues = []
        
        # Issue 1: Reasoning traces
        reasoning_patterns = [
            (r'<think>.*?</think>', "Think blocks"),
            (r'^Let me.*?(?=\n|\.|:|$)', "Let me reasoning"),
            (r'^I need to.*?(?=\n|\.|:|$)', "I need to reasoning"),
            (r'^First.*?(?=\n|\.|:|$)', "First reasoning"),
            (r'^Okay.*?(?=\n|\.|:|$)', "Okay reasoning"),
            (r'^Now.*?(?=\n|\.|:|$)', "Now reasoning"),
        ]
        
        for pattern, desc in reasoning_patterns:
            matches = re.finditer(pattern, text, re.DOTALL | re.MULTILINE)
            for match in matches:
                line_num = text[:match.start()].count('\n') + 1
                issues.append(LaTeXIssue(
                    issue_type="reasoning_trace",
                    description=f"{desc} found in output",
                    line_number=line_num,
                    context=match.group()[:50] + "..." if len(match.group()) > 50 else match.group(),
                    severity="medium"
                ))
        
        return issues
