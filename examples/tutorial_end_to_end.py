#!/usr/bin/env python3
"""Compatibility entry point for the installed, packaged synthetic walkthrough."""

from scrfu.tutorial import main, run_tutorial

__all__ = ["run_tutorial"]

if __name__ == "__main__":
    main()
