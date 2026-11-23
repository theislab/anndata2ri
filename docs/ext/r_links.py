"""Sphinx extension for links to R documentation."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, ClassVar
from urllib.error import HTTPError
from urllib.parse import urljoin
from urllib.request import urlopen

from docutils import nodes
from lxml import html
from sphinx.roles import XRefRole
from sphinx.util import logging


if TYPE_CHECKING:
    from sphinx.application import Sphinx
    from sphinx.environment import BuildEnvironment


logger = logging.getLogger(__name__)


class RManRefRole(XRefRole):
    """R reference role."""

    nodeclass: ClassVar[nodes.Node] = nodes.reference
    topic_cache: ClassVar[dict[str, dict[str, str]]] = {}
    """pkg → alias → url"""

    cls: bool

    def __init__(self, *a, cls: bool = False, **kw) -> None:  # noqa: ANN002, ANN003
        """Set self.cls."""
        super().__init__(*a, **kw)
        self.cls = cls

    def _get_man(self, pkg: str, alias: str) -> str:
        pkg_cache = type(self).topic_cache.setdefault(pkg)
        if not pkg_cache:
            for repo in ['R-patched', 'cran', 'bioc']:
                try:
                    pkg_cache = self._fetch_cache(repo, pkg)
                    break
                except HTTPError:
                    pass
            else:
                return None
            type(self).topic_cache[pkg] = pkg_cache
        return pkg_cache.get(alias)

    def _fetch_cache(self, repo: str, pkg: str) -> dict[str, str]:
        if not repo.startswith('R'):
            # R-universe builds current manuals for all CRAN and Bioconductor packages
            base = f'https://{repo}.r-universe.dev'
            with urlopen(f'{base}/api/packages/{pkg}') as con:  # noqa: S310
                help_pages = json.load(con)['_help']
            return {
                topic: f'{base}/{pkg}/doc/manual.html#{page["page"]}' for page in help_pages for topic in page['topics']
            }

        url = f'https://stat.ethz.ch/R-manual/{repo}/library/{pkg}/html/00Index.html'
        with urlopen(url) as con:
            txt = con.read().decode(con.headers.get_content_charset())
        doc = html.fromstring(txt)
        return {tr[0][0].text: urljoin(url, tr[0][0].attrib['href']) for tr in doc.xpath('//tr')}

    def process_link(
        self,
        env: BuildEnvironment,  # noqa: ARG002
        refnode: nodes.reference,
        has_explicit_title: bool,  # noqa: ARG002, FBT001
        title: str,
        target: str,
    ) -> tuple[str, str]:
        """Derive link title and URL from target."""
        qualified = not target.startswith('~')
        if not qualified:
            target = target[1:]
        package, symbol = target.split('::')
        title = target if qualified else symbol
        topic = symbol
        if self.cls:
            topic += '-class'
        url = self._get_man(package, topic)
        refnode['refuri'] = url
        if not url:
            logger.warning('R topic %s not found.', target, location=refnode)
        return title, url


def setup(app: Sphinx) -> None:
    """Set Sphinx extension up."""
    app.add_role('rman', RManRefRole())
    app.add_role('rcls', RManRefRole(cls=True))
