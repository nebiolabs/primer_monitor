import 'init_jquery';

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
export function registerPageModule(guardFn, initFn, cleanupFn) {
    let initialized = false;

    function maybeInit() {
        if (!guardFn()) return;
        if (initialized) return;
        initialized = true;
        initFn();
    }

    $(document).on('turbo:before-cache', function() {
        if (!guardFn()) return;
        initialized = false;
        cleanupFn?.();
    });

    $(document).on('turbo:load', maybeInit);

    // Handle first Turbo navigation to this page: turbo:load already fired before
    // this module loaded, so initialize now if the DOM is already ready.
    maybeInit();
}
