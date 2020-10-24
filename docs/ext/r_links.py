import logging
from typing import Tuple

from docutils import nodes
from sphinx.application import Sphinx
from sphinx.environment import BuildEnvironment
from sphinx.roles import XRefRole


class RManRefRole(XRefRole):
    nodeclass = nodes.reference

    topic_cache = {}

    def __init__(self, *a, cls: bool = False, **kw):
        super().__init__(*a, **kw)
        self.cls = cls

    def _get_man(self, pkg: str, alias: str):
        from urllib.error import HTTPError

        pkg_cache = type(self).topic_cache.setdefault(pkg)
        if not pkg_cache:
            for repo in ["R-patched", "cran", "bioc"]:
                try:
                    pkg_cache = self._fetch_cache(repo, pkg)
                    break
                except HTTPError:
                    pass
            else:
                return None
            type(self).topic_cache[pkg] = pkg_cache
        return pkg_cache.get(alias)

    def _fetch_cache(self, repo: str, pkg: str):
        from lxml import html
        from urllib.request import urlopen
        from urllib.parse import urljoin

        if repo.startswith("R"):
            url = f"https://stat.ethz.ch/R-manual/{repo}/library/{pkg}/html/00Index.html"
            tr_xpath = "//tr"
            get = lambda tr: (tr[0][0].text, tr[0][0].attrib["href"])
        else:
            url = f"https://rdrr.io/{repo}/{pkg}/api/"
            tr_xpath = "//div[@id='body-content']//tr[./td]"
            get = lambda tr: (tr[0].text, tr[1][0].attrib["href"])

        with urlopen(url) as con:
            txt = con.read().decode(con.headers.get_content_charset())
        doc = html.fromstring(txt)
        cache = {}
        for tr in doc.xpath(tr_xpath):
            topic, href = get(tr)
            cache[topic] = urljoin(url, href)
        return cache

    def process_link(
        self, env: BuildEnvironment, refnode: nodes.reference, has_explicit_title: bool, title: str, target: str
    ) -> Tuple[str, str]:
        qualified = not target.startswith("~")
        if not qualified:
            target = target[1:]
        package, symbol = target.split("::")
        title = target if qualified else symbol
        topic = symbol
        if self.cls:
            topic += "-class"
        url = self._get_man(package, topic)
        refnode["refuri"] = url
        if not url:
            logging.warning(f"R topic {target} not found.")
        return title, url

    # def result_nodes(self, document: nodes.document, env: BuildEnvironment, node: nodes.reference, is_ref: bool):
    #    target = node.get('reftarget')
    #    if target:
    #        node.attributes['refuri'] = target
    #    return [node], []


def setup(app: Sphinx):
    app.add_role("rman", RManRefRole())
    app.add_role("rcls", RManRefRole(cls=True))
