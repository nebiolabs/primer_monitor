// Registers a page-specific module that works correctly with Turbo navigation.
//
// The race condition: when a page-specific <script type="module"> is injected by
// Turbo during navigation, turbo:load fires before the module has loaded and
// registered its listener. guardFn/initFn/cleanupFn address this by also calling
// init immediately at module-load time.
//
// guardFn  - returns truthy when the current page owns this module
// initFn   - called once per page visit to initialize
// cleanupFn - called in turbo:before-cache to tear down before Turbo snapshots the page
//
// Initialization is tracked per <body>, because Turbo puts in a new body on every render. A form that redirects back
// to the page it was on (approving a primer set, for one) renders a new body without turbo:before-cache, so a plain
// initialized flag would skip setting up the new page.
export function registerPageModule(guardFn, initFn, cleanupFn) {
    let initializedBody = null;

    function maybeInit() {
        if (!guardFn()) return;
        if (initializedBody === document.body) return;
        // the previous page's setup was never torn down
        if (initializedBody) cleanupFn?.();
        initializedBody = document.body;
        initFn();
    }

    document.addEventListener('turbo:before-cache', function() {
        if (!guardFn() || initializedBody !== document.body) return;
        initializedBody = null;
        cleanupFn?.();
    });

    document.addEventListener('turbo:load', maybeInit);

    // Handle first Turbo navigation to this page: turbo:load already fired before
    // this module loaded, so initialize now if the DOM is already ready.
    maybeInit();
}
