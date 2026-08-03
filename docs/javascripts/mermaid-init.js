document$.subscribe(function() {
  if (typeof mermaid !== "undefined") {
    mermaid.initialize({
      startOnLoad: true,
      theme: document.body.getAttribute("data-md-color-scheme") === "slate" ? "dark" : "default"
    });
    mermaid.init(undefined, document.querySelectorAll(".mermaid"));
  }
});
