"""
LLM interface for the agent brain.

Supports multiple providers:
  - Gemini (Google, free tier) — default
  - Claude (Anthropic)
  - ChatGPT (OpenAI)

Usage:
  llm = create_llm(provider="gemini")            # uses env var
  llm = create_llm(provider="claude", api_key="sk-ant-...")
  llm = create_llm(provider="chatgpt", api_key="sk-...")
"""

import re
import json
import logging
from typing import Optional

from .config import (
    GEMINI_MODEL, CLAUDE_MODEL, CHATGPT_MODEL,
    get_gemini_api_key, get_anthropic_api_key, get_openai_api_key,
)


class BaseLLM:
    """Base class for LLM providers."""

    provider_name: str = "base"
    model_name: str = ""

    def __init__(self):
        self.available = False

    def ask(self, prompt: str) -> Optional[str]:
        raise NotImplementedError

    def ask_json(self, prompt: str) -> Optional[dict]:
        """Send prompt expecting JSON response. Parses and returns dict."""
        response = self.ask(prompt)
        if not response:
            return None
        try:
            clean = response.strip()
            clean = re.sub(r'^```json\s*', '', clean)
            clean = re.sub(r'\s*```$', '', clean)
            return json.loads(clean)
        except (json.JSONDecodeError, ValueError) as e:
            logging.warning(f"  Could not parse LLM JSON: {e}")
            logging.debug(f"  Raw: {response[:300]}")
            return None


class GeminiLLM(BaseLLM):
    """Google Gemini API (free tier).

    Install: pip install google-genai
    """

    provider_name = "gemini"
    model_name = GEMINI_MODEL

    def __init__(self, api_key: Optional[str] = None):
        super().__init__()
        self.client = None

        api_key = api_key or get_gemini_api_key()
        if not api_key:
            logging.info("LLM: Disabled (no GEMINI_API_KEY)")
            return

        try:
            from google import genai
            logging.getLogger("google.genai").setLevel(logging.WARNING)
            self.client = genai.Client(api_key=api_key)
            self.available = True
            logging.info(f"LLM: {self.model_name} (Gemini)")
        except ImportError:
            logging.warning("  google-genai not installed — run: pip install google-genai")
        except Exception as e:
            logging.warning(f"  Gemini init failed: {e}")

    def ask(self, prompt: str) -> Optional[str]:
        if not self.available:
            return None
        try:
            response = self.client.models.generate_content(
                model=self.model_name,
                contents=prompt,
            )
            return response.text
        except Exception as e:
            logging.warning(f"  Gemini API error: {e}")
            return None


class ClaudeLLM(BaseLLM):
    """Anthropic Claude API.

    Install: pip install anthropic
    """

    provider_name = "claude"
    model_name = CLAUDE_MODEL

    def __init__(self, api_key: Optional[str] = None):
        super().__init__()
        self.client = None

        api_key = api_key or get_anthropic_api_key()
        if not api_key:
            logging.info("LLM: Disabled (no ANTHROPIC_API_KEY)")
            return

        try:
            import anthropic
            self.client = anthropic.Anthropic(api_key=api_key)
            self.available = True
            logging.info(f"LLM: {self.model_name} (Claude)")
        except ImportError:
            logging.warning("  anthropic not installed — run: pip install anthropic")
        except Exception as e:
            logging.warning(f"  Claude init failed: {e}")

    def ask(self, prompt: str) -> Optional[str]:
        if not self.available:
            return None
        try:
            message = self.client.messages.create(
                model=self.model_name,
                max_tokens=4096,
                messages=[{"role": "user", "content": prompt}],
            )
            return message.content[0].text
        except Exception as e:
            logging.warning(f"  Claude API error: {e}")
            return None


class ChatGPTLLM(BaseLLM):
    """OpenAI ChatGPT API.

    Install: pip install openai
    """

    provider_name = "chatgpt"
    model_name = CHATGPT_MODEL

    def __init__(self, api_key: Optional[str] = None):
        super().__init__()
        self.client = None

        api_key = api_key or get_openai_api_key()
        if not api_key:
            logging.info("LLM: Disabled (no OPENAI_API_KEY)")
            return

        try:
            import openai
            self.client = openai.OpenAI(api_key=api_key)
            self.available = True
            logging.info(f"LLM: {self.model_name} (ChatGPT)")
        except ImportError:
            logging.warning("  openai not installed — run: pip install openai")
        except Exception as e:
            logging.warning(f"  ChatGPT init failed: {e}")

    def ask(self, prompt: str) -> Optional[str]:
        if not self.available:
            return None
        try:
            response = self.client.chat.completions.create(
                model=self.model_name,
                messages=[{"role": "user", "content": prompt}],
            )
            return response.choices[0].message.content
        except Exception as e:
            logging.warning(f"  ChatGPT API error: {e}")
            return None


_LLM_PROVIDERS = {
    "gemini": GeminiLLM,
    "claude": ClaudeLLM,
    "chatgpt": ChatGPTLLM,
}


def create_llm(provider: str = "gemini", api_key: Optional[str] = None) -> BaseLLM:
    """Factory to create an LLM instance for the given provider.

    Args:
        provider: One of "gemini", "claude", "chatgpt".
        api_key: Optional API key (overrides environment variable).

    Returns:
        An LLM instance (may have available=False if not configured).
    """
    cls = _LLM_PROVIDERS.get(provider)
    if cls is None:
        logging.warning(f"  Unknown LLM provider '{provider}', falling back to Gemini")
        cls = GeminiLLM
    return cls(api_key=api_key)
