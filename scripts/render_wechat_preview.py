#!/usr/bin/env python3
"""Compatibility entry point for the 16S-style WeChat bundle builder.

Use ``scripts/build_wechat_review_bundle.py`` in new workflows.  Keeping this
entry point prevents older local commands from invoking the retired blue/orange
header-card and hard-compaction renderer.
"""

from build_wechat_review_bundle import main


if __name__ == "__main__":
    raise SystemExit(main())
