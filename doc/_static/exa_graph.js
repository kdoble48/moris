/* Relational graph of the EXA examples (doc/_ext/exa_graph.py writes the data; D3 draws it).
   Each <div class="exa-graph" data-mode="ego|map"> holds its nodes and edges as JSON. */
(function () {
  "use strict";

  const CATEGORY_COLORS = {
    optimization: "#c8541c", structure: "#087f99", thermal: "#b8860b",
    fluid: "#4b7bec", mesh_generation: "#7a8a99",
  };
  const CATEGORY_LABELS = {
    optimization: "optimization", structure: "structure", thermal: "thermal",
    fluid: "fluid", mesh_generation: "mesh generation",
  };

  function color(node) {
    return CATEGORY_COLORS[node.category] || "#999";
  }

  function tooltipHtml(node, edge) {
    const tags = node.tags.map((t) => `<span class="exa-tag">${t}</span>`).join(" ");
    const shared = edge
      ? `<ul class="exa-graph-shared">${(edge.shared.length ? edge.shared : ["same category"]).map((s) => `<li>${s}</li>`).join("")}</ul>`
      : "";
    const thumb = node.thumb ? `<img src="${node.thumb}" alt="">` : "";
    return `${thumb}<strong>${node.title}</strong><div>${tags}</div>${shared}`;
  }

  function draw(container) {
    if (typeof d3 === "undefined") { return; }
    const data = JSON.parse(container.querySelector("script").textContent);
    const mode = container.dataset.mode;
    const width = container.clientWidth || 600;
    const height = container.clientHeight || 320;
    const centre = data.nodes.find((n) => n.centre);
    const byId = new Map(data.nodes.map((n) => [n.id, n]));
    const edgeOf = new Map(data.edges.map((e) => [`${e.source}|${e.target}`, e]));

    const svg = d3.select(container).append("svg")
      .attr("viewBox", `0 0 ${width} ${height}`).attr("width", "100%").attr("height", "100%");
    const layer = svg.append("g");
    const tip = d3.select(container).append("div").attr("class", "exa-graph-tip").style("display", "none");
    function placeTip(ev) {
      // keep the tooltip inside the box (overflow is hidden): flip it above / to the left near an edge
      const box = container.getBoundingClientRect();
      const w = tip.node().offsetWidth;
      const h = tip.node().offsetHeight;
      const x = ev.offsetX + 14 + w > box.width ? Math.max(0, ev.offsetX - 14 - w) : ev.offsetX + 14;
      const y = ev.offsetY + 14 + h > box.height ? Math.max(0, ev.offsetY - 14 - h) : ev.offsetY + 14;
      tip.style("left", `${x}px`).style("top", `${y}px`);
    }

    let zoom = null;
    if (mode === "map") {
      zoom = d3.zoom().scaleExtent([0.3, 4]).on("zoom", (ev) => layer.attr("transform", ev.transform));
      svg.call(zoom);
      const bar = d3.select(container).append("div").attr("class", "exa-graph-bar");
      const tags = Array.from(new Set(data.nodes.flatMap((n) => n.tags))).sort();
      const select = bar.append("select");
      select.append("option").attr("value", "").text("highlight a feature…");
      tags.forEach((t) => select.append("option").attr("value", t).text(t));
      select.on("change", function () { tagFilter = this.value; applyFilters(); });
      const legend = bar.append("span").attr("class", "exa-graph-legend");
      Object.keys(CATEGORY_LABELS).forEach((c) => {
        legend.append("button").attr("type", "button").attr("class", "exa-graph-family").attr("data-category", c)
          .attr("title", "click to hide or show this family")
          .html(`<i style="background:${CATEGORY_COLORS[c]}"></i>${CATEGORY_LABELS[c]}`)
          .on("click", function () {
            if (hidden.has(c)) { hidden.delete(c); } else { hidden.add(c); }
            d3.select(this).classed("exa-graph-off", hidden.has(c));
            applyFilters();
          });
      });
    }
    const hidden = new Set();
    let tagFilter = "";
    function applyFilters() {
      node.style("display", (n) => (hidden.has(n.category) ? "none" : null))
        .classed("exa-graph-dim", (n) => tagFilter && !n.tags.includes(tagFilter));
      const shown = (e) => !hidden.has(e.source.category) && !hidden.has(e.target.category);
      link.style("display", (e) => (shown(e) ? null : "none"))
        .classed("exa-graph-dim", (e) => tagFilter && !(e.source.tags.includes(tagFilter) && e.target.tags.includes(tagFilter)));
      linkHit.style("display", (e) => (shown(e) ? null : "none"));
    }

    const simulation = d3.forceSimulation(data.nodes)
      .force("link", d3.forceLink(data.edges).id((n) => n.id)
        .distance((e) => (mode === "ego" ? 170 : 70) * (1.5 - Math.min(e.weight, 1)))
        .strength((e) => 0.3 + 0.7 * Math.min(e.weight, 1)))
      .force("charge", d3.forceManyBody().strength(mode === "ego" ? -900 : -220))
      .force("collide", d3.forceCollide(mode === "ego" ? 60 : 22))
      .force("center", d3.forceCenter(width / 2, height / 2))
      .force("x", d3.forceX(width / 2).strength(mode === "ego" ? 0.04 : 0.12))
      .force("y", d3.forceY(height / 2).strength(mode === "ego" ? 0.06 : 0.14));
    if (centre) { centre.fx = width / 2; centre.fy = height / 2; }

    const link = layer.append("g").selectAll("line").data(data.edges).join("line")
      .attr("class", "exa-graph-link")
      .attr("stroke-width", (e) => (mode === "ego" ? 1 + 5 : 0.5 + 2.5) * Math.min(e.weight, 1));

    const linkHit = layer.append("g").selectAll("line").data(data.edges).join("line")
      .attr("class", "exa-graph-link-hit")
      .on("mousemove", (ev, e) => {
        const shared = e.shared.length ? e.shared.map((s) => `<li>${s}</li>`).join("") : "<li>same category</li>";
        tip.style("display", "block")
          .html(`<strong>${e.source.title}</strong> &harr; <strong>${e.target.title}</strong>`
            + `<ul class="exa-graph-shared">${shared}</ul><div class="exa-graph-score">similarity ${e.weight.toFixed(2)}</div>`)
;
        placeTip(ev);
      })
      .on("mouseleave", () => tip.style("display", "none"));

    const node = layer.append("g").selectAll("g").data(data.nodes).join("g")
      .attr("class", (n) => "exa-graph-node" + (n.centre ? " exa-graph-centre" : ""))
      .style("cursor", (n) => (n.url && !n.centre ? "pointer" : "default"))
      .on("click", (ev, n) => { if (n.url && !n.centre) { window.location.href = n.url; } })
      .on("mousemove", (ev, n) => {
        const edge = centre && !n.centre
          ? (edgeOf.get(`${centre.id}|${n.id}`) || edgeOf.get(`${n.id}|${centre.id}`)) : null;
        tip.style("display", "block").html(tooltipHtml(n, edge))
;
        placeTip(ev);
      })
      .on("mouseleave", () => tip.style("display", "none"))
      .call(d3.drag()
        .on("start", (ev, n) => { if (!ev.active) simulation.alphaTarget(0.3).restart(); n.fx = n.x; n.fy = n.y; })
        .on("drag", (ev, n) => { n.fx = ev.x; n.fy = ev.y; })
        .on("end", (ev, n) => { if (!ev.active) simulation.alphaTarget(0); if (!n.centre) { n.fx = null; n.fy = null; } }));

    node.append("circle").attr("r", (n) => (n.centre ? 14 : mode === "ego" ? 10 : 7)).attr("fill", color);
    const label = node.append("text").text((n) => n.title)
      .attr("class", mode === "map" ? "exa-graph-label exa-graph-label-small" : "exa-graph-label");
    function placeLabels() {
      label.attr("dx", (n) => (mode === "ego" && !n.centre ? (n.x < width / 2 ? -14 : 14) : 0))
        .attr("dy", (n) => (mode === "ego" && !n.centre ? 4 : n.centre ? 30 : 18))
        .attr("text-anchor", (n) => (mode === "ego" && !n.centre ? (n.x < width / 2 ? "end" : "start") : "middle"));
    }

    simulation.stop();
    for (let i = 0; i < 300; i += 1) { simulation.tick(); }
    // fit the settled layout to the box, with room for the labels
    const padX = mode === "ego" ? 150 : 90;
    const padY = mode === "ego" ? 40 : 50;
    const xs = data.nodes.map((n) => n.x);
    const ys = data.nodes.map((n) => n.y);
    const spanX = Math.max(1, Math.max(...xs) - Math.min(...xs));
    const spanY = Math.max(1, Math.max(...ys) - Math.min(...ys));
    const scale = Math.min((width - 2 * padX) / spanX, (height - 2 * padY) / spanY, 1.4);
    const fit = d3.zoomIdentity
      .translate(width / 2 - scale * (Math.min(...xs) + spanX / 2), height / 2 - scale * (Math.min(...ys) + spanY / 2))
      .scale(scale);
    if (zoom) { svg.call(zoom.transform, fit); } else { layer.attr("transform", fit); }
    simulation.on("tick", () => {
      link.attr("x1", (e) => e.source.x).attr("y1", (e) => e.source.y)
        .attr("x2", (e) => e.target.x).attr("y2", (e) => e.target.y);
      linkHit.attr("x1", (e) => e.source.x).attr("y1", (e) => e.source.y)
        .attr("x2", (e) => e.target.x).attr("y2", (e) => e.target.y);
      node.attr("transform", (n) => `translate(${n.x},${n.y})`);
      placeLabels();
    });
    simulation.tick();
    simulation.on("tick")();
  }

  function init() {
    document.querySelectorAll(".exa-graph").forEach(draw);
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
