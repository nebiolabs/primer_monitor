import igv from "igv";

// The page's one igv.js browser, created and removed through igv's API so Turbo navigation never leaves a stale
// or duplicate browser behind. igv.createBrowser is async, so removal waits for a browser that is still loading.
let pending = null;

export function createBrowser(div, config) {
    removeBrowser();
    pending = igv.createBrowser(div, config);
    return pending;
}

export function removeBrowser() {
    if (!pending) return;
    const loading = pending;
    pending = null;
    loading.then(browser => igv.removeBrowser(browser)).catch(() => {});
}
