"""Sphinx extension for links to R documentation."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, ClassVar
from urllib.error import HTTPError
from urllib.parse import urljoin
from urllib.request import urlopen

from docutils import nodes
from lxml import html
from sphinx.roles import XRefRole


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
        if repo.startswith('R'):
            url = f'https://stat.ethz.ch/R-manual/{repo}/library/{pkg}/html/00Index.html'
            tr_xpath = '//tr'

            def get(tr: html.HtmlElement) -> tuple[str, str]:
                return tr[0][0].text, tr[0][0].attrib['href']

        else:
            url = f'https://rdrr.io/{repo}/{pkg}/api/'
            tr_xpath = "//div[@id='body-content']//tr[./td]"

            def get(tr: html.HtmlElement) -> tuple[str, str]:
                return tr[0].text, tr[1][0].attrib['href']

        with urlopen(url) as con:  # noqa: S310
            txt = con.read().decode(con.headers.get_content_charset())
        doc = html.fromstring(txt)
        cache = {}
        for tr in doc.xpath(tr_xpath):
            topic, href = get(tr)
            cache[topic] = urljoin(url, href)
        return cache

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
            logger.warning('R topic %s not found.', target)
        return title, url


def setup(app: Sphinx) -> None:
    """Set Sphinx extension up."""
    app.add_role('rman', RManRefRole())
    app.add_role('rcls', RManRefRole(cls=True))
