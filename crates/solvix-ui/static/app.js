/* global document, window, navigator, ResizeObserver */
"use strict";

(function () {
  const embeddedModel = {
    ok: true,
    job: "cantilever",
    nodes: [
      [1, 0.0, 0.0],
      [2, 0.0, 0.05],
    ],
    elements: [[1, "cps4", [1, 2, 5, 4]]],
    nsets: { FIX: [1, 2, 3], TIP: [41, 82, 123, 164, 205] },
    elsets: { ALL: [1, 2, 3] },
    materials: [{ name: "STEEL", young: 210000000000.0, poisson: 0.3 }],
    bcs: [{ nodes: [1, 2, 3], dofs: [1, 2] }],
    cloads: [[205, 2, -200.0]],
  };

  const embeddedResults = {
    ok: true,
    job: "cantilever",
    displacements: [
      [0.0, 0.0],
      [0.000004, -0.000008],
    ],
    von_mises: [120000000.0, 119000000.0],
    max_disp_mm: 13.832,
    max_von_mises_mpa: 79.483,
    iterations: 390,
    residual: 6.933e-11,
  };

  const embeddedModel3D = {
    ok: true,
    job: "cantilever3d",
    nodes: [
      [1,0.0,0.0,0.0], [2,1.0,0.0,0.0], [3,2.0,0.0,0.0], [4,3.0,0.0,0.0], [5,4.0,0.0,0.0],
      [6,0.0,1.0,0.0], [7,1.0,1.0,0.0], [8,2.0,1.0,0.0], [9,3.0,1.0,0.0], [10,4.0,1.0,0.0],
      [11,0.0,0.0,0.5], [12,1.0,0.0,0.5], [13,2.0,0.0,0.5], [14,3.0,0.0,0.5], [15,4.0,0.0,0.5],
      [16,0.0,1.0,0.5], [17,1.0,1.0,0.5], [18,2.0,1.0,0.5], [19,3.0,1.0,0.5], [20,4.0,1.0,0.5],
      [21,0.0,0.0,1.0], [22,1.0,0.0,1.0], [23,2.0,0.0,1.0], [24,3.0,0.0,1.0], [25,4.0,0.0,1.0],
      [26,0.0,1.0,1.0], [27,1.0,1.0,1.0], [28,2.0,1.0,1.0], [29,3.0,1.0,1.0], [30,4.0,1.0,1.0],
    ],
    elements: [
      [1,"c3d8",[1,2,7,6,11,12,17,16]],
      [2,"c3d8",[2,3,8,7,12,13,18,17]],
      [3,"c3d8",[3,4,9,8,13,14,19,18]],
      [4,"c3d8",[4,5,10,9,14,15,20,19]],
      [5,"c3d8",[11,12,17,16,21,22,27,26]],
      [6,"c3d8",[12,13,18,17,22,23,28,27]],
      [7,"c3d8",[13,14,19,18,23,24,29,28]],
      [8,"c3d8",[14,15,20,19,24,25,30,29]],
    ],
    nsets: {
      FIX: [1,6,11,16,21,26],
      TIP: [5,10,15,20,25,30],
    },
    elsets: { ALL: [1,2,3,4,5,6,7,8] },
    materials: [{"name": "STEEL","young": 210000000000.0,"poisson": 0.3}],
    bcs: [{"nodes": [1,6,11,16,21,26],"dofs": [1,2,3]}],
    cloads: [
      [5,2,-1666.6666666666667], [10,2,-1666.6666666666667], [15,2,-1666.6666666666667],
      [20,2,-1666.6666666666667], [25,2,-1666.6666666666667], [30,2,-1666.6666666666667]
    ],
  };

  const embeddedResults3D = {
    ok: true,
    job: "cantilever3d",
    displacements: [
      [-0.0,-0.0,0.0], [-1.0448979591836738e-11,-4.476190476190476e-06,0.0], [-4.179591836734695e-11,-8.38095238095238e-06,0.0],
      [-9.404081632653064e-11,-1.1142857142857143e-05,0.0], [-1.671836734693878e-10,-1.219047619047619e-05,0.0], [-0.0,-0.0,0.0],
      [-1.0448979591836738e-11,-4.476190476190476e-06,0.0], [-4.179591836734695e-11,-8.38095238095238e-06,0.0], [-9.404081632653064e-11,-1.1142857142857143e-05,0.0],
      [-1.671836734693878e-10,-1.219047619047619e-05,0.0], [-0.0,-0.0,0.0], [-1.0448979591836738e-11,-4.476190476190476e-06,0.0],
      [-4.179591836734695e-11,-8.38095238095238e-06,0.0], [-9.404081632653064e-11,-1.1142857142857143e-05,0.0], [-1.671836734693878e-10,-1.219047619047619e-05,0.0],
      [-0.0,-0.0,0.0], [-1.0448979591836738e-11,-4.476190476190476e-06,0.0], [-4.179591836734695e-11,-8.38095238095238e-06,0.0],
      [-9.404081632653064e-11,-1.1142857142857143e-05,0.0], [-1.671836734693878e-10,-1.219047619047619e-05,0.0], [-0.0,-0.0,0.0],
      [-1.0448979591836738e-11,-4.476190476190476e-06,0.0], [-4.179591836734695e-11,-8.38095238095238e-06,0.0], [-9.404081632653064e-11,-1.1142857142857143e-05,0.0],
      [-1.671836734693878e-10,-1.219047619047619e-05,0.0], [-0.0,-0.0,0.0], [-1.0448979591836738e-11,-4.476190476190476e-06,0.0],
      [-4.179591836734695e-11,-8.38095238095238e-06,0.0], [-9.404081632653064e-11,-1.1142857142857143e-05,0.0], [-1.671836734693878e-10,-1.219047619047619e-05,0.0]
    ],
    von_mises: [238000.0,171000.0,114000.0,68000.0,221000.0,155000.0,101000.0,57000.0],
    max_disp_mm: 0.01219047619047619,
    max_von_mises_mpa: 0.238,
    iterations: 96,
    residual: 4.8e-11,
  };

  const scriptUrl =
    (document.currentScript && document.currentScript.src) ||
    new URL("app.js", document.baseURI).href;
  const mockModelUrl = new URL("mock/model.json", scriptUrl);
  const mockResultsUrl = new URL("mock/results.json", scriptUrl);
  const mockModel3DUrl = new URL("mock/model3d.json", scriptUrl);
  const mockResults3DUrl = new URL("mock/results3d.json", scriptUrl);

  const dom = {};
  const requiredIds = [
    "data-source", "job-form", "deck-path", "run-button", "job-select",
    "node-count", "element-count", "element-types", "set-count", "set-list",
    "material-count", "material-list", "bc-count", "bc-list", "load-count",
    "load-list", "step-count", "step-list", "job-tree-count", "job-list",
    "show-node-dots", "show-node-ids", "show-element-ids", "highlight-sets",
    "show-bcs", "show-loads", "show-undeformed-overlay", "viewport",
    "viewport-mode-label",
    "fit-button", "wireframe-button", "contour-u-button",
    "contour-vm-button", "undeformed-overlay-button", "deform-slider",
    "deform-output", "result-chips", "console-output", "console-form",
    "command-input", "command-hint", "status-job", "status-mesh",
    "status-cursor", "status-zoom", "status-message",
  ];

  for (const id of requiredIds) {
    const found = document.getElementById(id);
    if (!found) {
      throw new Error(`Missing required element #${id}`);
    }
    dom[id] = found;
  }

  const canvas = dom.viewport;
  const ctx = canvas.getContext("2d", { alpha: false });

  const state = {
    module: "model",
    job: "cantilever",
    model: null,
    modelSource: "not loaded",
    modelError: null,
    results: null,
    resultsSource: "not loaded",
    resultsError: null,
    runResponse: null,
    jobs: ["cantilever3d", "cantilever"],
    is3D: false,
    modelBounds: null,
    camera: { scale: 100, cx: 0, cy: 0 },
    camera3D: {
      scale: 100,
      cx: 0,
      cy: 0,
      cz: 0,
      rotation: initial3DRotation(),
    },
    pointer: null,
    cursorWorld: null,
    display: {
      showNodeDots: false,
      showNodeIds: false,
      showElementIds: false,
      highlightSets: false,
      showBCs: true,
      showLoads: true,
      wireframe: true,
      undeformedOverlay: true,
      contourMode: "none",
      deformationPercent: 5,
    },
    selectedSets: new Map(),
    consoleEntries: [],
    commandHistory: [],
    historyIndex: -1,
    warnings: [],
  };

  function number(value, fallback = 0) {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : fallback;
  }

  function escapeHTML(value) {
    return String(value).replace(/[&<>'"]/g, (ch) => {
      return {
        "&": "&amp;",
        "<": "&lt;",
        ">": "&gt;",
        "'": "&#39;",
        "\"": "&quot;",
      }[ch];
    });
  }

  function pretty(value) {
    return JSON.stringify(value, null, 2);
  }

  function fmt(value, digits = 3) {
    if (!Number.isFinite(value)) {
      return "—";
    }
    const abs = Math.abs(value);
    if (abs !== 0 && (abs < 0.001 || abs >= 1000000)) {
      return value.toExponential(Math.max(2, digits - 1));
    }
    return value.toFixed(digits);
  }

  function fmtLength(meters) {
    const abs = Math.abs(meters);
    if (abs >= 1) return `${fmt(meters, 3)} m`;
    if (abs >= 0.001 || abs === 0) return `${fmt(meters * 1000, 3)} mm`;
    if (abs >= 0.000001) return `${fmt(meters * 1e6, 3)} µm`;
    return `${fmt(meters * 1e9, 3)} nm`;
  }

  const HEX_EDGE_INDICES = Object.freeze([
    [0, 1], [1, 2], [2, 3], [3, 0],
    [4, 5], [5, 6], [6, 7], [7, 4],
    [0, 4], [1, 5], [2, 6], [3, 7],
  ].map(Object.freeze));

  const TET_EDGE_INDICES = Object.freeze([
    [0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3],
  ].map(Object.freeze));

  function is3DModel(nodes, elements) {
    for (const node of nodes.values()) {
      if (Math.abs(node.z) > 1e-12) return true;
    }
    for (const element of elements.values()) {
      if (element.type === "c3d8" || element.type === "c3d4") return true;
    }
    return false;
  }

  function reconstructC3D8Depth(nodes, elements) {
    // Compatibility for API payloads emitted before issue #6: nodes are
    // [id,x,y] even for C3D8 decks. Abaqus hex connectivity gives a bottom
    // face (0..3) followed by its corresponding top face (4..7), so a stack's
    // missing depth coordinate can be recovered as a one-dimensional level
    // chain. The in-plane extent supplies a conservative layer spacing.
    if (nodes.size === 0 || elements.size === 0) return nodes;
    if (!Array.from(elements.values()).every((element) => element.type === "c3d8")) {
      return nodes;
    }

    const parent = new Map();
    const child = new Map();
    for (const element of elements.values()) {
      if (element.connectivity.length !== 8) return nodes;
      for (let i = 0; i < 4; i += 1) {
        const lower = element.connectivity[i];
        const upper = element.connectivity[i + 4];
        if (!nodes.has(lower) || !nodes.has(upper) || lower === upper) return nodes;
        if (parent.has(upper) && parent.get(upper) !== lower) return nodes;
        if (child.has(lower) && child.get(lower) !== upper) return nodes;
        parent.set(upper, lower);
        child.set(lower, upper);
      }
    }

    const levels = new Map();
    const visit = (nodeId, level) => {
      const known = levels.get(nodeId);
      if (known !== undefined) return known === level;
      levels.set(nodeId, level);
      const next = child.get(nodeId);
      return next === undefined || visit(next, level + 1);
    };

    const roots = Array.from(nodes.keys()).filter((nodeId) => !parent.has(nodeId));
    if (!roots.length || !roots.every((nodeId) => visit(nodeId, 0))) return nodes;
    if (levels.size !== nodes.size) return nodes;

    let maxLevel = 0;
    levels.forEach((level) => {
      maxLevel = Math.max(maxLevel, level);
    });
    if (maxLevel === 0) return nodes;

    const xs = Array.from(nodes.values(), (node) => node.x);
    const ys = Array.from(nodes.values(), (node) => node.y);
    const range = Math.min(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys));
    const first = elements.values().next().value;
    const edgeLengths = first.connectivity.slice(0, 4).map((id, index, ids) => {
      const a = nodes.get(id);
      const b = nodes.get(ids[(index + 1) % 4]);
      return Math.hypot(a.x - b.x, a.y - b.y);
    }).sort((a, b) => a - b);
    const fallbackSpacing = edgeLengths[Math.floor(edgeLengths.length / 2)] || 1;
    const spacing = range > 0 ? range / maxLevel : fallbackSpacing;
    const centerLevel = maxLevel / 2;

    const output = new Map();
    nodes.forEach((node, id) => {
      output.set(id, Object.assign({}, node, {
        z: (levels.get(id) - centerLevel) * spacing,
      }));
    });
    return output;
  }

  function normalizeModel(raw) {
    const nodes = new Map();
    const elements = new Map();
    const nsets = new Map();
    const elsets = new Map();
    let hasDepthCoordinate = false;

    for (const item of Array.isArray(raw.nodes) ? raw.nodes : []) {
      const id = Number(item && item[0]);
      if (!Number.isFinite(id)) continue;
      nodes.set(id, {
        id,
        x: number(item[1]),
        y: number(item[2]),
        z: number(item[3]),
      });
      hasDepthCoordinate = hasDepthCoordinate || item.length >= 4;
    }

    for (const item of Array.isArray(raw.elements) ? raw.elements : []) {
      const id = Number(item && item[0]);
      if (!Number.isFinite(id) || !Array.isArray(item[2])) continue;
      elements.set(id, {
        id,
        type: String(item[1] || "unknown").toLowerCase(),
        connectivity: item[2].map((nodeId) => Number(nodeId)),
      });
    }

    const physicalNodes = hasDepthCoordinate
      ? nodes
      : reconstructC3D8Depth(nodes, elements);

    const setMapper = (value) => new Set(Array.isArray(value) ? value.map(Number) : []);
    Object.entries(raw.nsets || {}).forEach(([name, value]) => nsets.set(name, setMapper(value)));
    Object.entries(raw.elsets || {}).forEach(([name, value]) => elsets.set(name, setMapper(value)));

    return {
      job: String(raw.job || state.job),
      nodes: physicalNodes,
      elements,
      nsets,
      elsets,
      materials: Array.isArray(raw.materials) ? raw.materials : [],
      bcs: Array.isArray(raw.bcs) ? raw.bcs : [],
      cloads: Array.isArray(raw.cloads) ? raw.cloads : [],
      steps: Array.isArray(raw.steps) ? raw.steps : [],
      jobs: Array.isArray(raw.jobs) ? raw.jobs : [],
      raw,
    };
  }

  function normalizeResults(raw, model) {
    const displacements = (Array.isArray(raw.displacements) ? raw.displacements : [])
      .map((item) => {
        if (Array.isArray(item)) {
          return [number(item[0]), number(item[1]), number(item[2])];
        }
        if (item && typeof item === "object") {
          return [number(item.x), number(item.y), number(item.z)];
        }
        return [0, 0, 0];
      });

    const vonMises = (Array.isArray(raw.von_mises) ? raw.von_mises : [])
      .map((item) => number(item));

    return {
      job: String(raw.job || (model ? model.job : state.job)),
      displacements,
      vonMises,
      raw,
    };
  }

  function boundsOf(nodes) {
    if (!nodes || nodes.size === 0) return null;
    const values = Array.from(nodes.values());
    const bounds = values.reduce((acc, node) => ({
      minX: Math.min(acc.minX, node.x),
      maxX: Math.max(acc.maxX, node.x),
      minY: Math.min(acc.minY, node.y),
      maxY: Math.max(acc.maxY, node.y),
      minZ: Math.min(acc.minZ, node.z),
      maxZ: Math.max(acc.maxZ, node.z),
    }), {
      minX: Infinity, maxX: -Infinity, minY: Infinity, maxY: -Infinity,
      minZ: Infinity, maxZ: -Infinity,
    });
    return bounds;
  }

  function unionBounds(a, b) {
    if (!a) return b;
    if (!b) return a;
    return {
      minX: Math.min(a.minX, b.minX),
      maxX: Math.max(a.maxX, b.maxX),
      minY: Math.min(a.minY, b.minY),
      maxY: Math.max(a.maxY, b.maxY),
      minZ: Math.min(a.minZ, b.minZ),
      maxZ: Math.max(a.maxZ, b.maxZ),
    };
  }

  function displayNodes() {
    if (!state.model) return new Map();
    const output = new Map();
    const factor = deformationFactor();
    const displacement = state.results ? state.results.displacements : [];
    let index = 0;
    for (const node of state.model.nodes.values()) {
      const d = displacement[index] || [0, 0, 0];
      output.set(node.id, {
        id: node.id,
        x: node.x + d[0] * factor,
        y: node.y + d[1] * factor,
        z: node.z + d[2] * factor,
      });
      index += 1;
    }
    return output;
  }

  function maxDisplacementMagnitude() {
    if (!state.results) return 0;
    return state.results.displacements.reduce((maximum, item) => {
      return Math.max(maximum, Math.hypot(item[0], item[1], item[2] || 0));
    }, 0);
  }

  function deformationFactor() {
    const bounds = state.modelBounds;
    if (!bounds || !state.results || state.display.deformationPercent === 0) return 0;
    const maximum = maxDisplacementMagnitude();
    if (!(maximum > 0)) return 0;
    const diagonal = state.is3D
      ? Math.hypot(bounds.maxX - bounds.minX, bounds.maxY - bounds.minY, bounds.maxZ - bounds.minZ)
      : Math.hypot(bounds.maxX - bounds.minX, bounds.maxY - bounds.minY);
    if (!(diagonal > 0)) return 0;
    return (state.display.deformationPercent / 100) * diagonal / maximum;
  }

  function matrixIdentity() {
    return [1, 0, 0, 0, 1, 0, 0, 0, 1];
  }

  function matrixMultiply(a, b) {
    return [
      a[0] * b[0] + a[1] * b[3] + a[2] * b[6],
      a[0] * b[1] + a[1] * b[4] + a[2] * b[7],
      a[0] * b[2] + a[1] * b[5] + a[2] * b[8],
      a[3] * b[0] + a[4] * b[3] + a[5] * b[6],
      a[3] * b[1] + a[4] * b[4] + a[5] * b[7],
      a[3] * b[2] + a[4] * b[5] + a[5] * b[8],
      a[6] * b[0] + a[7] * b[3] + a[8] * b[6],
      a[6] * b[1] + a[7] * b[4] + a[8] * b[7],
      a[6] * b[2] + a[7] * b[5] + a[8] * b[8],
    ];
  }

  function matrixTranspose(matrix) {
    return [
      matrix[0], matrix[3], matrix[6],
      matrix[1], matrix[4], matrix[7],
      matrix[2], matrix[5], matrix[8],
    ];
  }

  function matrixVectorMultiply(matrix, vector) {
    return [
      matrix[0] * vector[0] + matrix[1] * vector[1] + matrix[2] * vector[2],
      matrix[3] * vector[0] + matrix[4] * vector[1] + matrix[5] * vector[2],
      matrix[6] * vector[0] + matrix[7] * vector[1] + matrix[8] * vector[2],
    ];
  }

  function rotationAxisAngle(axis, angle) {
    const length = Math.hypot(axis[0], axis[1], axis[2]);
    if (!(length > 0) || angle === 0) return matrixIdentity();
    const [x, y, z] = axis.map((value) => value / length);
    const cosine = Math.cos(angle);
    const sine = Math.sin(angle);
    const oneMinusCosine = 1 - cosine;
    return [
      cosine + x * x * oneMinusCosine,
      x * y * oneMinusCosine - z * sine,
      x * z * oneMinusCosine + y * sine,
      y * x * oneMinusCosine + z * sine,
      cosine + y * y * oneMinusCosine,
      y * z * oneMinusCosine - x * sine,
      z * x * oneMinusCosine - y * sine,
      z * y * oneMinusCosine + x * sine,
      cosine + z * z * oneMinusCosine,
    ];
  }

  function rotationBetweenUnitVectors(a, b) {
    const dot = Math.min(1, Math.max(-1, a[0] * b[0] + a[1] * b[1] + a[2] * b[2]));
    const axis = [
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
    ];
    const axisLength = Math.hypot(...axis);
    if (axisLength > 1e-8) return rotationAxisAngle(axis, Math.atan2(axisLength, dot));
    if (dot > 0) return matrixIdentity();
    // A 180-degree rotation is axis-ambiguous. Choose any perpendicular axis.
    const seed = Math.abs(a[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    const perpendicular = [
      a[1] * seed[2] - a[2] * seed[1],
      a[2] * seed[0] - a[0] * seed[2],
      a[0] * seed[1] - a[1] * seed[0],
    ];
    return rotationAxisAngle(perpendicular, Math.PI);
  }

  function rotationXY(xAngle, yAngle) {
    return matrixMultiply(
      rotationAxisAngle([0, 1, 0], yAngle),
      rotationAxisAngle([1, 0, 0], xAngle)
    );
  }

  function initial3DRotation() {
    return rotationXY(-22 * Math.PI / 180, -32 * Math.PI / 180);
  }

  function arcballVector(x, y, width, height) {
    const radius = Math.max(1, Math.min(width, height) / 2);
    const vx = (2 * x - width) / radius;
    const vy = (height - 2 * y) / radius;
    const lengthSquared = vx * vx + vy * vy;
    return [vx, vy, lengthSquared <= 1 ? Math.sqrt(1 - lengthSquared) : 0];
  }

  function projectOrthographic(node, camera, rotation, width, height) {
    const delta = [node.x - camera.cx, node.y - camera.cy, node.z - camera.cz];
    const view = matrixVectorMultiply(rotation, delta);
    return {
      x: width / 2 + view[0] * camera.scale,
      y: height / 2 - view[1] * camera.scale,
      depth: view[2],
    };
  }

  function fitOrthographic(bounds, rotation, width, height, margin = 72, legendWidth = 82) {
    if (!bounds) {
      return { scale: 100, cx: 0, cy: 0, cz: 0, rotation };
    }
    const center = {
      x: (bounds.minX + bounds.maxX) / 2,
      y: (bounds.minY + bounds.maxY) / 2,
      z: (bounds.minZ + bounds.maxZ) / 2,
    };
    const camera = { scale: 1, cx: center.x, cy: center.y, cz: center.z, rotation };
    let minX = Infinity;
    let maxX = -Infinity;
    let minY = Infinity;
    let maxY = -Infinity;
    const xUnit = [rotation[0], rotation[3], rotation[6]];
    const yUnit = [rotation[1], rotation[4], rotation[7]];
    const corner = (dx, dy, dz) => {
      const point = { x: center.x + dx, y: center.y + dy, z: center.z + dz };
      const px = (point.x - center.x) * xUnit[0]
        + (point.y - center.y) * xUnit[1]
        + (point.z - center.z) * xUnit[2];
      const py = (point.x - center.x) * yUnit[0]
        + (point.y - center.y) * yUnit[1]
        + (point.z - center.z) * yUnit[2];
      minX = Math.min(minX, px);
      maxX = Math.max(maxX, px);
      minY = Math.min(minY, py);
      maxY = Math.max(maxY, py);
    };
    const halfX = Math.max((bounds.maxX - bounds.minX) / 2, 0);
    const halfY = Math.max((bounds.maxY - bounds.minY) / 2, 0);
    const halfZ = Math.max((bounds.maxZ - bounds.minZ) / 2, 0);
    [-1, 1].forEach((sx) => [-1, 1].forEach((sy) => [-1, 1].forEach((sz) => {
      corner(sx * halfX, sy * halfY, sz * halfZ);
    })));
    const availableWidth = Math.max(30, width - margin * 2 - legendWidth);
    const availableHeight = Math.max(30, height - margin * 2);
    const scale = Math.min(
      availableWidth / Math.max(maxX - minX, 1e-12),
      availableHeight / Math.max(maxY - minY, 1e-12)
    );
    return Object.assign({}, camera, { scale: Math.min(1e10, Math.max(1e-6, scale)) });
  }

  function screenToCenterWorld(x, y, camera, rotation, width, height) {
    const cameraPoint = [
      (x - width / 2) / camera.scale,
      (height / 2 - y) / camera.scale,
      0,
    ];
    const worldDelta = matrixVectorMultiply(matrixTranspose(rotation), cameraPoint);
    return {
      x: camera.cx + worldDelta[0],
      y: camera.cy + worldDelta[1],
      z: camera.cz + worldDelta[2],
    };
  }

  function panOrthographic(camera, rotation, dxScreen, dyScreen) {
    const cameraDelta = [-dxScreen / camera.scale, dyScreen / camera.scale, 0];
    const worldDelta = matrixVectorMultiply(matrixTranspose(rotation), cameraDelta);
    return {
      scale: camera.scale,
      cx: camera.cx + worldDelta[0],
      cy: camera.cy + worldDelta[1],
      cz: camera.cz + worldDelta[2],
      rotation,
    };
  }

  function zoomOrthographicAt(camera, rotation, x, y, factor, width, height) {
    const scale = Math.min(1e10, Math.max(1e-6, camera.scale * factor));
    const target = screenToCenterWorld(x, y, camera, rotation, width, height);
    const nextCamera = { scale, cx: target.x, cy: target.y, cz: target.z, rotation };
    const after = screenToCenterWorld(x, y, nextCamera, rotation, width, height);
    return {
      scale,
      cx: nextCamera.cx + target.x - after.x,
      cy: nextCamera.cy + target.y - after.y,
      cz: nextCamera.cz + target.z - after.z,
      rotation,
    };
  }

  function edgeIndicesForElement(type) {
    if (type === "c3d8") return HEX_EDGE_INDICES;
    if (type === "c3d4") return TET_EDGE_INDICES;
    return [];
  }

  function build3DEdgeRecords(elements, projectedNodes, elementValues = []) {
    const records = [];
    let elementIndex = 0;
    for (const element of elements.values()) {
      const indices = edgeIndicesForElement(element.type);
      for (const [a, b] of indices) {
        const start = projectedNodes.get(element.connectivity[a]);
        const end = projectedNodes.get(element.connectivity[b]);
        if (!start || !end) continue;
        records.push({
          elementId: element.id,
          value: elementValues[elementIndex],
          start,
          end,
          depth: (start.depth + end.depth) / 2,
        });
      }
      elementIndex += 1;
    }
    return records;
  }

  function painterSortEdges(records) {
    return records.slice().sort((a, b) => a.depth - b.depth);
  }

  async function readJSON(response) {
    const text = await response.text();
    let payload;
    try {
      payload = JSON.parse(text);
    } catch (_) {
      throw new Error(`HTTP ${response.status}: non-JSON response`);
    }
    if (!response.ok) {
      const detail = payload && payload.error ? `: ${payload.error}` : "";
      throw new Error(`HTTP ${response.status}${detail}`);
    }
    return payload;
  }

  async function requestJSON(path, init) {
    const response = await fetch(path, init);
    return readJSON(response);
  }

  async function loadMock(url, embedded) {
    try {
      return await requestJSON(url);
    } catch (error) {
      if (location.protocol === "file:") {
        return embedded;
      }
      throw error;
    }
  }

  function is3DJob(job) {
    return /3d$/i.test(String(job || ""));
  }

  function mockModelFor(job) {
    return is3DJob(job) ? mockModel3DUrl : mockModelUrl;
  }

  function mockResultsFor(job) {
    return is3DJob(job) ? mockResults3DUrl : mockResultsUrl;
  }

  function embeddedModelFor(job) {
    return is3DJob(job) ? embeddedModel3D : embeddedModel;
  }

  function embeddedResultsFor(job) {
    return is3DJob(job) ? embeddedResults3D : embeddedResults;
  }

  function transportError(error) {
    return {
      ok: false,
      stage: "transport",
      error: String(error && error.message ? error.message : error),
    };
  }

  async function loadModel(job, options = {}) {
    const requestedJob = String(job || "cantilever");
    let source = "api";
    try {
      const payload = await requestJSON(`/api/model?job=${encodeURIComponent(requestedJob)}`);
      if (payload.ok === false) {
        state.modelError = payload;
        source = "api-error";
        if (!options.silent) logJSON("error", payload, "GET model");
        return false;
      }
      applyModel(payload, source);
      if (!options.silent) logJSON("ok", payload, "GET model");
      return true;
    } catch (apiError) {
      try {
        const payload = await loadMock(mockModelFor(requestedJob), embeddedModelFor(requestedJob));
        applyModel(payload, "mock");
        state.modelError = transportError(apiError);
        if (!options.silent) {
          logJSON("info", state.modelError, "model API fallback");
          logJSON("ok", payload, "GET mock model");
        }
        return true;
      } catch (mockError) {
        const payload = transportError(mockError);
        state.modelError = payload;
        if (!options.silent) logJSON("error", payload, "GET model");
        setStatus("Unable to load model");
        return false;
      }
    }
  }

  async function loadResults(job, options = {}) {
    const requestedJob = String(job || state.job || "cantilever");
    let source = "api";
    try {
      const payload = await requestJSON(`/api/results/${encodeURIComponent(requestedJob)}`);
      if (payload.ok === false) {
        state.resultsError = payload;
        state.results = null;
        if (!options.silent) logJSON("error", payload, "GET results");
        updateInterface();
        render();
        return false;
      }
      applyResults(payload, source);
      if (!options.silent) logJSON("ok", payload, "GET results");
      return true;
    } catch (apiError) {
      // Never splice a compact 3D mock result onto a full API-sized model.
      if (state.modelSource === "api") {
        state.results = null;
        state.resultsError = transportError(apiError);
        updateInterface();
        render();
        return false;
      }
      try {
        const payload = await loadMock(mockResultsFor(requestedJob), embeddedResultsFor(requestedJob));
        applyResults(payload, "mock");
        state.resultsError = transportError(apiError);
        if (!options.silent) {
          logJSON("info", state.resultsError, "results API fallback");
          logJSON("ok", payload, "GET mock results");
        }
        return true;
      } catch (mockError) {
        const payload = transportError(mockError);
        state.results = null;
        state.resultsError = payload;
        if (!options.silent) logJSON("error", payload, "GET results");
        updateInterface();
        render();
        return false;
      }
    }
  }

  function applyModel(raw, source) {
    state.model = normalizeModel(raw);
    state.modelSource = source;
    state.modelError = null;
    state.job = state.model.job;
    state.is3D = is3DModel(state.model.nodes, state.model.elements);
    state.modelBounds = boundsOf(state.model.nodes);
    state.selectedSets.clear();
    state.camera3D.rotation = initial3DRotation();
    canvas.setAttribute("aria-label", state.is3D
      ? "Three-dimensional finite element model viewport"
      : "Finite element model viewport");
    if (!state.jobs.includes(state.job)) {
      state.jobs.push(state.job);
    }
    ensureJobOption(state.job);
    fitView();
    updateInterface();
    render();
    setStatus(`Loaded model ${state.job}`);
  }

  function applyResults(raw, source) {
    state.results = normalizeResults(raw, state.model);
    state.resultsSource = source;
    state.resultsError = null;
    if (!state.jobs.includes(state.results.job)) {
      state.jobs.push(state.results.job);
    }
    ensureJobOption(state.results.job);
    state.job = state.results.job;
    state.warnings = [];
    if (state.model && state.results.displacements.length !== state.model.nodes.size) {
      state.warnings.push(
        `Result has ${state.results.displacements.length} displacement rows; model has ${state.model.nodes.size} nodes.`
      );
    }
    if (state.model && state.results.vonMises.length !== state.model.elements.size) {
      state.warnings.push(
        `Result has ${state.results.vonMises.length} von Mises values; model has ${state.model.elements.size} elements.`
      );
    }
    updateInterface();
    render();
  }

  function ensureJobOption(job) {
    let option = Array.from(dom["job-select"].options).find((item) => item.value === job);
    if (!option) {
      option = document.createElement("option");
      option.value = job;
      option.textContent = job;
      dom["job-select"].appendChild(option);
    }
    dom["job-select"].value = job;
  }

  function updateInterface() {
    const model = state.model;
    dom["node-count"].textContent = model ? model.nodes.size : 0;
    dom["element-count"].textContent = model ? model.elements.size : 0;
    dom["set-count"].textContent = model ? model.nsets.size + model.elsets.size : 0;
    dom["material-count"].textContent = model ? model.materials.length : 0;
    dom["bc-count"].textContent = model ? model.bcs.length : 0;
    dom["load-count"].textContent = model ? model.cloads.length : 0;
    dom["step-count"].textContent = model ? model.steps.length : 0;
    dom["job-tree-count"].textContent = state.jobs.length;
    dom["status-job"].textContent = `job: ${state.job}`;
    const types = model
      ? Array.from(model.elements.values()).reduce((counts, element) => {
        counts[element.type] = (counts[element.type] || 0) + 1;
        return counts;
      }, {})
      : {};
    const typeText = Object.keys(types).length
      ? Object.entries(types).map(([type, count]) => `${type.toUpperCase()} × ${count}`).join(", ")
      : "no element kinds";
    dom["status-mesh"].textContent = model
      ? `nodes: ${model.nodes.size} · elements: ${model.elements.size} · ${
        state.is3D ? "3D" : "2D"
      } · ${typeText}`
      : "nodes: 0 · elements: 0";
    dom["data-source"].textContent = `model: ${state.modelSource} · results: ${state.resultsSource}`;
    dom["viewport-mode-label"].textContent = viewportLabel();

    dom["element-types"].textContent = Object.keys(types).length
      ? typeText
      : "No elements";

    renderSetList();
    renderSimpleLists();
    renderResultChips();
    updateToolbarStates();
  }

  function renderSetList() {
    const model = state.model;
    dom["set-list"].replaceChildren();
    if (!model || (model.nsets.size === 0 && model.elsets.size === 0)) {
      const item = document.createElement("div");
      item.className = "muted-line";
      item.style.paddingLeft = "21px";
      item.textContent = "No sets in model payload";
      dom["set-list"].appendChild(item);
      return;
    }

    const appendSet = (kind, name, set) => {
      const row = document.createElement("div");
      row.className = "set-row";
      const label = document.createElement("label");
      const input = document.createElement("input");
      input.type = "checkbox";
      input.checked = state.selectedSets.get(name) === kind;
      input.addEventListener("change", () => {
        if (input.checked) {
          state.selectedSets.set(name, kind);
          state.display.highlightSets = true;
          dom["highlight-sets"].checked = true;
        } else {
          state.selectedSets.delete(name);
        }
        render();
        updateInterface();
      });
      const nameSpan = document.createElement("span");
      nameSpan.textContent = name;
      const kindSpan = document.createElement("span");
      kindSpan.className = "set-kind";
      kindSpan.textContent = `${kind} · ${set.size}`;
      label.append(input, nameSpan, kindSpan);
      row.appendChild(label);
      dom["set-list"].appendChild(row);
    };

    model.nsets.forEach((set, name) => appendSet("nset", name, set));
    model.elsets.forEach((set, name) => appendSet("elset", name, set));
  }

  function renderSimpleLists() {
    const model = state.model;
    const materialList = model ? model.materials : [];
    dom["material-list"].replaceChildren(...materialList.map((material) => {
      const li = document.createElement("li");
      const name = material && material.name ? material.name : "unnamed";
      const young = number(material && material.young);
      const poisson = number(material && material.poisson);
      li.textContent = `${name} · E=${fmt(young / 1e9, 1)} GPa · ν=${fmt(poisson, 2)}`;
      li.title = pretty(material);
      return li;
    }));
    if (!materialList.length) {
      dom["material-list"].replaceChildren(emptyItem("No materials"));
    }

    const bcList = model ? model.bcs : [];
    dom["bc-list"].replaceChildren(...bcList.map((bc) => {
      const li = document.createElement("li");
      const nodes = Array.isArray(bc && bc.nodes) ? bc.nodes : [];
      const dofs = Array.isArray(bc && bc.dofs) ? bc.dofs : [];
      li.textContent = `${nodes.length} nodes · DOF ${dofs.join(", ") || "—"}`;
      li.title = pretty(bc);
      return li;
    }));
    if (!bcList.length) dom["bc-list"].replaceChildren(emptyItem("No BCs"));

    const loadList = model ? model.cloads : [];
    dom["load-list"].replaceChildren(...loadList.map((load) => {
      const li = document.createElement("li");
      if (Array.isArray(load)) {
        li.textContent = `node ${load[0]} · DOF ${load[1]} · ${fmt(number(load[2]), 2)} N`;
      } else {
        li.textContent = "CLOAD";
      }
      li.title = pretty(load);
      return li;
    }));
    if (!loadList.length) dom["load-list"].replaceChildren(emptyItem("No loads"));

    const stepList = model ? model.steps : [];
    dom["step-list"].replaceChildren(...stepList.map((step) => {
      const li = document.createElement("li");
      li.textContent = typeof step === "string" ? step : pretty(step);
      li.title = pretty(step);
      return li;
    }));
    if (!stepList.length) dom["step-list"].replaceChildren(emptyItem("No step data"));

    dom["job-list"].replaceChildren(...state.jobs.map((job) => {
      const li = document.createElement("li");
      const button = document.createElement("button");
      button.type = "button";
      button.className = "mini-copy";
      button.textContent = job;
      button.addEventListener("click", () => {
        dom["job-select"].value = job;
        switchJob(job);
      });
      li.appendChild(button);
      return li;
    }));
  }

  function emptyItem(text) {
    const li = document.createElement("li");
    li.className = "muted-line";
    li.textContent = text;
    return li;
  }

  function renderResultChips() {
    const chips = [];
    if (state.results) {
      const raw = state.results.raw;
      const maxVM = raw && Number.isFinite(Number(raw.max_von_mises_mpa))
        ? Number(raw.max_von_mises_mpa)
        : Math.max(0, ...state.results.vonMises) / 1e6;
      const maxDisp = raw && Number.isFinite(Number(raw.max_disp_mm))
        ? Number(raw.max_disp_mm)
        : maxDisplacementMagnitude() * 1000;
      chips.push(`max |U| ${fmt(maxDisp, 3)} mm`);
      chips.push(`max VM ${fmt(maxVM, 3)} MPa`);
      if (Number.isFinite(Number(raw && raw.iterations))) {
        chips.push(`${Number(raw.iterations)} iterations`);
      }
    } else {
      chips.push("No result");
    }
    if (state.display.contourMode !== "none") {
      chips.push(`contour ${state.display.contourMode === "vm" ? "von Mises" : "|U|"}`);
    }
    chips.push(`deform ${state.display.deformationPercent.toFixed(1)}%`);
    dom["result-chips"].replaceChildren(...chips.map((text) => {
      const span = document.createElement("span");
      span.className = "chip";
      span.textContent = text;
      return span;
    }));
  }

  function updateToolbarStates() {
    dom["wireframe-button"].setAttribute("aria-pressed", String(state.display.wireframe));
    dom["contour-u-button"].setAttribute("aria-pressed", String(state.display.contourMode === "u"));
    dom["contour-vm-button"].setAttribute("aria-pressed", String(state.display.contourMode === "vm"));
    dom["undeformed-overlay-button"].setAttribute("aria-pressed", String(state.display.undeformedOverlay));
    dom["deform-slider"].value = String(state.display.deformationPercent);
    dom["deform-output"].value = `${state.display.deformationPercent.toFixed(1)}% bbox`;
    document.querySelectorAll(".module-button").forEach((button) => {
      button.classList.toggle("is-active", button.dataset.module === state.module);
    });
  }

  function viewportLabel() {
    const projection = state.is3D ? "3D · " : "";
    if (!state.results) return `${projection}Undeformed`;
    const contour = state.display.contourMode === "vm"
      ? "von Mises"
      : state.display.contourMode === "u" ? "|U|" : "no contour";
    return `${projection}Deformed · ${contour}`;
  }

  function setStatus(message) {
    dom["status-message"].textContent = message;
  }

  function logJSON(kind, payload, label) {
    const entry = document.createElement("article");
    entry.className = "console-entry";
    entry.dataset.kind = kind;
    const meta = document.createElement("div");
    meta.className = "console-meta";
    const strong = document.createElement("strong");
    strong.textContent = label || kind.toUpperCase();
    const time = document.createElement("time");
    time.dateTime = new Date().toISOString();
    time.textContent = new Date().toLocaleTimeString();
    meta.append(strong, time);

    const pre = document.createElement("pre");
    pre.className = "console-json";
    pre.textContent = pretty(payload);

    const copy = document.createElement("button");
    copy.type = "button";
    copy.className = "mini-copy entry-copy";
    copy.textContent = "Copy";
    copy.addEventListener("click", () => copyToClipboard(payload, label || "console entry"));

    entry.append(meta, pre, copy);
    dom["console-output"].appendChild(entry);
    dom["console-output"].scrollTop = dom["console-output"].scrollHeight;
    state.consoleEntries.push({ kind, label: label || kind, payload, at: time.dateTime });
  }

  async function copyToClipboard(payload, label) {
    const text = pretty(payload);
    const legacyCopy = async () => {
      const textarea = document.createElement("textarea");
      textarea.value = text;
      textarea.style.position = "fixed";
      textarea.style.opacity = "0";
      document.body.appendChild(textarea);
      textarea.select();
      const copied = document.execCommand("copy");
      textarea.remove();
      if (!copied) throw new Error("clipboard command returned false");
    };
    try {
      if (navigator.clipboard && window.isSecureContext) {
        await navigator.clipboard.writeText(text);
      } else {
        await legacyCopy();
      }
      setStatus(`Copied ${label} JSON`);
    } catch (error) {
      try {
        await legacyCopy();
        setStatus(`Copied ${label} JSON`);
      } catch (fallbackError) {
        setStatus(`Copy failed: ${fallbackError.message}`);
      }
    }
  }

  function payloadForCopy(kind) {
    const model = state.model;
    switch (kind) {
      case "model":
        return model ? model.raw : { ok: false, error: "model not loaded" };
      case "results":
        return state.results
          ? state.results.raw
          : { ok: false, error: "results not loaded" };
      case "tree.nodes":
        return model ? { nodes: model.raw.nodes || [] } : { nodes: [] };
      case "tree.elements":
        return model ? { elements: model.raw.elements || [] } : { elements: [] };
      case "tree.sets":
        return model
          ? { nsets: model.raw.nsets || {}, elsets: model.raw.elsets || {} }
          : { nsets: {}, elsets: {} };
      case "tree.materials":
        return model ? { materials: model.raw.materials || [] } : { materials: [] };
      case "tree.bcs":
        return model ? { bcs: model.raw.bcs || [] } : { bcs: [] };
      case "tree.loads":
        return model ? { cloads: model.raw.cloads || [] } : { cloads: [] };
      case "tree.steps":
        return model ? { steps: model.raw.steps || [] } : { steps: [] };
      case "tree.jobs":
        return { jobs: state.jobs, active: state.job };
      case "job":
        return {
          job: state.job,
          jobs: state.jobs,
          deckPath: dom["deck-path"].value,
          lastRun: state.runResponse,
          modelSource: state.modelSource,
          resultsSource: state.resultsSource,
        };
      case "viewport":
        return {
          job: state.job,
          mode: viewportLabel(),
          module: state.module,
          display: Object.assign({}, state.display, {
            selectedSets: Object.fromEntries(state.selectedSets),
          }),
          renderer: state.is3D ? "3d-orthographic" : "2d",
          camera: Object.assign({}, state.camera),
          camera3D: Object.assign({}, state.camera3D),
          canvas: { width: canvas.clientWidth, height: canvas.clientHeight },
          warnings: state.warnings,
          modelSource: state.modelSource,
          resultsSource: state.resultsSource,
        };
      case "console":
        return { entries: state.consoleEntries };
      default:
        return { ok: false, error: `unknown copy source: ${kind}` };
    }
  }

  document.addEventListener("click", (event) => {
    const button = event.target.closest("[data-copy]");
    if (!button) return;
    event.preventDefault();
    const kind = button.dataset.copy;
    copyToClipboard(payloadForCopy(kind), kind);
  });

  /* Canvas geometry, interaction and drawing */
  function canvasSize() {
    return { width: Math.max(1, canvas.clientWidth), height: Math.max(1, canvas.clientHeight) };
  }

  function resizeCanvas() {
    const dpr = window.devicePixelRatio || 1;
    const { width, height } = canvasSize();
    canvas.width = Math.round(width * dpr);
    canvas.height = Math.round(height * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    render();
  }

  function worldToScreen(node) {
    const { width, height } = canvasSize();
    return {
      x: (node.x - state.camera.cx) * state.camera.scale + width / 2,
      y: (state.camera.cy - node.y) * state.camera.scale + height / 2,
    };
  }

  function screenToWorld(x, y) {
    const { width, height } = canvasSize();
    return {
      x: (x - width / 2) / state.camera.scale + state.camera.cx,
      y: state.camera.cy - (y - height / 2) / state.camera.scale,
    };
  }

  function fitView() {
    const bounds = currentFitBounds();
    const { width, height } = canvasSize();
    if (state.is3D) {
      state.camera3D = fitOrthographic(
        bounds,
        state.camera3D.rotation,
        width,
        height,
        62,
        state.display.contourMode === "none" ? 0 : 82
      );
      updateZoomStatus();
      return;
    }
    if (!bounds) {
      state.camera = { scale: 100, cx: 0, cy: 0 };
      updateZoomStatus();
      return;
    }
    const dx = Math.max(bounds.maxX - bounds.minX, 1e-12);
    const dy = Math.max(bounds.maxY - bounds.minY, 1e-12);
    const margin = 62;
    const legendMargin = state.display.contourMode === "none" ? 0 : 82;
    const availableWidth = Math.max(30, width - margin * 2 - legendMargin);
    const availableHeight = Math.max(30, height - margin * 2);
    state.camera = {
      scale: Math.min(availableWidth / dx, availableHeight / dy),
      cx: (bounds.minX + bounds.maxX) / 2,
      cy: (bounds.minY + bounds.maxY) / 2,
    };
    updateZoomStatus();
  }

  function currentFitBounds() {
    if (!state.model) return null;
    if (!state.results || state.display.deformationPercent === 0) {
      return state.modelBounds;
    }
    return unionBounds(state.modelBounds, boundsOf(displayNodes()));
  }

  function updateZoomStatus() {
    const scale = state.is3D ? state.camera3D.scale : state.camera.scale;
    dom["status-zoom"].textContent = `zoom: ${fmt(scale, 1)} px/m`;
  }

  canvas.addEventListener("pointerdown", (event) => {
    canvas.setPointerCapture(event.pointerId);
    const rect = canvas.getBoundingClientRect();
    const pan = event.shiftKey || event.button === 1 || event.button === 2;
    const localX = event.clientX - rect.left;
    const localY = event.clientY - rect.top;
    const { width, height } = canvasSize();
    state.pointer = {
      id: event.pointerId,
      x: event.clientX,
      y: event.clientY,
      mode: state.is3D && !pan ? "rotate" : "pan",
      vector: arcballVector(localX, localY, width, height),
      camera: Object.assign({}, state.is3D ? state.camera3D : state.camera),
      rotation: state.camera3D.rotation.slice(),
    };
    event.preventDefault();
    canvas.classList.add(state.pointer.mode === "rotate" ? "is-rotating" : "is-panning");
  });

  canvas.addEventListener("pointermove", (event) => {
    const rect = canvas.getBoundingClientRect();
    const localX = event.clientX - rect.left;
    const localY = event.clientY - rect.top;
    const { width, height } = canvasSize();
    if (state.is3D) {
      const world = screenToCenterWorld(
        localX,
        localY,
        state.camera3D,
        state.camera3D.rotation,
        width,
        height
      );
      state.cursorWorld = world;
      dom["status-cursor"].textContent =
        `x: ${fmt(world.x, 5)} y: ${fmt(world.y, 5)} z: ${fmt(world.z, 5)}`;
    } else {
      const world = screenToWorld(localX, localY);
      state.cursorWorld = world;
      dom["status-cursor"].textContent = `x: ${fmt(world.x, 5)} y: ${fmt(world.y, 5)}`;
    }

    if (!state.pointer || state.pointer.id !== event.pointerId) return;
    event.preventDefault();
    const dxScreen = event.clientX - state.pointer.x;
    const dyScreen = event.clientY - state.pointer.y;
    if (state.pointer.mode === "rotate") {
      const nextVector = arcballVector(localX, localY, width, height);
      const delta = rotationBetweenUnitVectors(state.pointer.vector, nextVector);
      state.camera3D.rotation = matrixMultiply(delta, state.pointer.rotation);
    } else if (state.is3D) {
      state.camera3D = panOrthographic(
        state.pointer.camera,
        state.pointer.rotation,
        dxScreen,
        dyScreen
      );
    } else {
      state.camera.cx = state.pointer.camera.cx - dxScreen / state.camera.scale;
      state.camera.cy = state.pointer.camera.cy + dyScreen / state.camera.scale;
    }
    render();
  });

  function finishPointer(event) {
    if (!state.pointer || state.pointer.id !== event.pointerId) return;
    if (canvas.hasPointerCapture(event.pointerId)) {
      canvas.releasePointerCapture(event.pointerId);
    }
    state.pointer = null;
    canvas.classList.remove("is-panning", "is-rotating");
  }

  canvas.addEventListener("pointerup", finishPointer);
  canvas.addEventListener("pointercancel", finishPointer);
  canvas.addEventListener("pointerleave", () => {
    state.cursorWorld = null;
    dom["status-cursor"].textContent = state.is3D ? "x: — y: — z: —" : "x: — y: —";
  });

  canvas.addEventListener("contextmenu", (event) => event.preventDefault());

  canvas.addEventListener("wheel", (event) => {
    event.preventDefault();
    const rect = canvas.getBoundingClientRect();
    const px = event.clientX - rect.left;
    const py = event.clientY - rect.top;
    const factor = Math.exp(-event.deltaY * 0.0014);
    const { width, height } = canvasSize();
    if (state.is3D) {
      state.camera3D = zoomOrthographicAt(
        state.camera3D,
        state.camera3D.rotation,
        px,
        py,
        factor,
        width,
        height
      );
    } else {
      const nextScale = Math.min(1e10, Math.max(1e-6, state.camera.scale * factor));
      const before = screenToWorld(px, py);
      state.camera.scale = nextScale;
      // Keep the world point under the cursor fixed as the scale changes.
      state.camera.cx = before.x - (px - width / 2) / nextScale;
      state.camera.cy = before.y + (py - height / 2) / nextScale;
    }
    updateZoomStatus();
    render();
  }, { passive: false });

  function jetColor(t) {
    const value = Math.min(1, Math.max(0, t));
    const stops = [
      [0.0, [0.0, 0.0, 0.5]],
      [0.12, [0.0, 0.0, 1.0]],
      [0.37, [0.0, 1.0, 1.0]],
      [0.62, [0.5, 1.0, 0.0]],
      [0.87, [1.0, 0.35, 0.0]],
      [1.0, [1.0, 0.0, 0.0]],
    ];
    for (let i = 1; i < stops.length; i += 1) {
      if (value <= stops[i][0]) {
        const [t0, c0] = stops[i - 1];
        const [t1, c1] = stops[i];
        const local = (value - t0) / (t1 - t0);
        return `rgb(${Math.round(255 * (c0[0] + (c1[0] - c0[0]) * local))}, ${
          Math.round(255 * (c0[1] + (c1[1] - c0[1]) * local))
        }, ${Math.round(255 * (c0[2] + (c1[2] - c0[2]) * local))})`;
      }
    }
    return "rgb(255,0,0)";
  }

  function contourDomain() {
    if (!state.results || state.display.contourMode === "none") return null;
    if (state.display.contourMode === "vm") {
      const values = state.results.vonMises;
      if (!values.length) return null;
      return { min: Math.min(...values), max: Math.max(...values), mode: "vm" };
    }
    const values = state.results.displacements.map((item) => Math.hypot(item[0], item[1], item[2] || 0));
    if (!values.length) return null;
    return { min: Math.min(...values), max: Math.max(...values), mode: "u" };
  }

  function normalizedColor(value, domain) {
    if (!domain) return jetColor(0);
    const range = domain.max - domain.min;
    return jetColor(range > 0 ? (value - domain.min) / range : 0.5);
  }

  function polygonFromConnectivity(connectivity, nodes) {
    const points = [];
    for (const nodeId of connectivity) {
      const node = nodes.get(Number(nodeId));
      if (node) points.push(node);
    }
    return points;
  }

  function projected3DNodes(nodes) {
    const { width, height } = canvasSize();
    const output = new Map();
    nodes.forEach((node, id) => {
      output.set(id, Object.assign(
        { id },
        projectOrthographic(node, state.camera3D, state.camera3D.rotation, width, height)
      ));
    });
    return output;
  }

  function draw3DMesh(nodes, deformed) {
    if (!state.model) return new Map();
    const projected = projected3DNodes(nodes);
    const domain = contourDomain();
    const elementValues = domain && domain.mode === "vm" ? state.results.vonMises : [];
    const records = painterSortEdges(
      build3DEdgeRecords(state.model.elements, projected, elementValues)
    );

    ctx.lineJoin = "round";
    ctx.lineCap = "round";
    ctx.lineWidth = 1.05;
    for (const edge of records) {
      if (domain && domain.mode === "vm" && Number.isFinite(edge.value)) {
        ctx.strokeStyle = normalizedColor(edge.value, domain);
        ctx.globalAlpha = 0.94;
      } else {
        ctx.strokeStyle = deformed ? "rgba(228, 232, 246, 0.9)" : "#d8c671";
        ctx.globalAlpha = deformed ? 0.88 : 0.94;
      }
      ctx.beginPath();
      ctx.moveTo(edge.start.x, edge.start.y);
      ctx.lineTo(edge.end.x, edge.end.y);
      ctx.stroke();
    }
    ctx.globalAlpha = 1;

    if (state.display.showElementIds) {
      ctx.fillStyle = "rgba(224, 229, 246, 0.84)";
      ctx.font = "9px ui-monospace, monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      for (const element of state.model.elements.values()) {
        let x = 0;
        let y = 0;
        let count = 0;
        for (const nodeId of element.connectivity) {
          const point = projected.get(nodeId);
          if (!point) continue;
          x += point.x;
          y += point.y;
          count += 1;
        }
        if (count) ctx.fillText(String(element.id), x / count, y / count);
      }
    }
    return projected;
  }

  function draw3DUndeformedOverlay() {
    if (!state.model || !state.results || !state.display.undeformedOverlay) return;
    if (deformationFactor() === 0) return;
    const projected = projected3DNodes(state.model.nodes);
    const records = painterSortEdges(build3DEdgeRecords(state.model.elements, projected));
    ctx.save();
    ctx.strokeStyle = "rgba(216, 198, 113, 0.58)";
    ctx.setLineDash([4, 3]);
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (const edge of records) {
      ctx.moveTo(edge.start.x, edge.start.y);
      ctx.lineTo(edge.end.x, edge.end.y);
    }
    ctx.stroke();
    ctx.restore();
  }

  function draw3DSelectedSets(projected) {
    if (!state.model || !state.display.highlightSets || state.selectedSets.size === 0) return;
    ctx.save();
    for (const [name, kind] of state.selectedSets) {
      const set = kind === "nset" ? state.model.nsets.get(name) : state.model.elsets.get(name);
      if (!set) continue;
      ctx.fillStyle = "rgba(0, 211, 255, 0.18)";
      ctx.strokeStyle = "rgba(0, 211, 255, 0.9)";
      ctx.lineWidth = 1.7;
      if (kind === "nset") {
        for (const nodeId of set) {
          const point = projected.get(nodeId);
          if (!point) continue;
          ctx.beginPath();
          ctx.arc(point.x, point.y, 6, 0, Math.PI * 2);
          ctx.fill();
          ctx.stroke();
        }
      } else {
        const records = build3DEdgeRecords(
          new Map(Array.from(state.model.elements).filter(([id]) => set.has(id))),
          projected
        );
        ctx.beginPath();
        for (const edge of records) {
          ctx.moveTo(edge.start.x, edge.start.y);
          ctx.lineTo(edge.end.x, edge.end.y);
        }
        ctx.stroke();
      }
    }
    ctx.restore();
  }

  function draw3DNodesAndIds(projected) {
    ctx.textAlign = "left";
    ctx.textBaseline = "alphabetic";
    for (const [id, point] of projected) {
      if (state.display.showNodeDots) {
        ctx.fillStyle = "#f2f4fb";
        ctx.beginPath();
        ctx.arc(point.x, point.y, 2.1, 0, Math.PI * 2);
        ctx.fill();
      }
      if (state.display.showNodeIds) {
        ctx.fillStyle = "rgba(184, 203, 255, 0.92)";
        ctx.font = "9px ui-monospace, monospace";
        ctx.fillText(String(id), point.x + 5, point.y - 5);
      }
    }
  }

  function draw3DNodalContour(projected) {
    const domain = contourDomain();
    if (!domain || domain.mode !== "u" || !state.results) return;
    const values = state.results.displacements.map((item) => Math.hypot(item[0], item[1], item[2] || 0));
    let index = 0;
    for (const point of Array.from(projected.values()).sort((a, b) => a.depth - b.depth)) {
      const value = values[index] || 0;
      ctx.fillStyle = normalizedColor(value, domain);
      ctx.beginPath();
      ctx.arc(point.x, point.y, 2.7, 0, Math.PI * 2);
      ctx.fill();
      index += 1;
    }
  }

  function draw3DBCsAndLoads(nodes, projected) {
    if (!state.model) return;
    if (state.display.showBCs) {
      ctx.fillStyle = "#ff922b";
      ctx.strokeStyle = "rgba(255, 146, 43, 0.45)";
      for (const bc of state.model.bcs) {
        const ids = Array.isArray(bc && bc.nodes) ? bc.nodes : [];
        for (const nodeId of ids) {
          const point = projected.get(Number(nodeId));
          if (!point) continue;
          ctx.fillRect(point.x - 4, point.y - 4, 8, 8);
          ctx.strokeRect(point.x - 4.5, point.y - 4.5, 9, 9);
        }
      }
    }

    if (!state.display.showLoads) return;
    const rotation = state.camera3D.rotation;
    ctx.strokeStyle = "#ff4747";
    ctx.fillStyle = "#ff4747";
    ctx.lineWidth = 2;
    for (const load of state.model.cloads) {
      if (!Array.isArray(load) || load.length < 3) continue;
      const point = projected.get(Number(load[0]));
      const node = nodes.get(Number(load[0]));
      if (!point || !node) continue;
      const dof = Number(load[1]);
      const axis = dof === 1 ? [1, 0, 0] : dof === 2 ? [0, 1, 0] : dof === 3 ? [0, 0, 1] : null;
      if (!axis) continue;
      const sign = Math.sign(number(load[2])) || 1;
      const viewDirection = matrixVectorMultiply(rotation, axis);
      const dx = viewDirection[0] * sign;
      const dy = -viewDirection[1] * sign;
      const length = Math.hypot(dx, dy);
      if (!(length > 1e-9)) continue;
      const end = {
        x: point.x + (dx / length) * 28,
        y: point.y + (dy / length) * 28,
      };
      const angle = Math.atan2(end.y - point.y, end.x - point.x);
      ctx.beginPath();
      ctx.moveTo(point.x, point.y);
      ctx.lineTo(end.x, end.y);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(end.x, end.y);
      ctx.lineTo(end.x - 7 * Math.cos(angle - 0.42), end.y - 7 * Math.sin(angle - 0.42));
      ctx.lineTo(end.x - 7 * Math.cos(angle + 0.42), end.y - 7 * Math.sin(angle + 0.42));
      ctx.closePath();
      ctx.fill();
    }
  }

  function triangleValue(p1, p2, p3, point) {
    const det = (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y);
    if (Math.abs(det) < 1e-18) return (p1.value + p2.value + p3.value) / 3;
    const lambda1 = ((p2.x - point.x) * (p3.y - point.y) - (p3.x - point.x) * (p2.y - point.y)) / det;
    const lambda2 = ((p3.x - point.x) * (p1.y - point.y) - (p1.x - point.x) * (p3.y - point.y)) / det;
    const lambda3 = 1 - lambda1 - lambda2;
    return p1.value * lambda1 + p2.value * lambda2 + p3.value * lambda3;
  }

  function drawMesh(nodes, deformed) {
    if (!state.model) return;
    const domain = contourDomain();
    const displacementValues = state.results
      ? state.results.displacements.map((item) => Math.hypot(item[0], item[1], item[2] || 0))
      : [];

    ctx.lineJoin = "round";
    ctx.lineCap = "round";

    let elementIndex = 0;
    const nodeOrder = new Map(
      Array.from(state.model.nodes.keys()).map((nodeId, index) => [nodeId, index])
    );
    for (const element of state.model.elements.values()) {
      const physical = polygonFromConnectivity(
        element.connectivity,
        deformed ? nodes : state.model.nodes
      );
      if (physical.length < 2) {
        elementIndex += 1;
        continue;
      }

      if (domain && physical.length >= 3) {
        if (domain.mode === "vm") {
          const value = state.results.vonMises[elementIndex];
          if (Number.isFinite(value)) {
            ctx.fillStyle = normalizedColor(value, domain);
            pathPolygon(physical);
            ctx.fill();
          }
        } else {
          drawNodalContour(physical, displacementValues, nodeOrder);
        }
      }

      const hasFilledContour = Boolean(domain && physical.length >= 3);
      if (state.display.wireframe || !hasFilledContour) {
        ctx.strokeStyle = deformed ? "rgba(226, 230, 244, 0.88)" : "#d8c671";
        ctx.lineWidth = deformed ? 1.0 : 1.15;
        pathPolygon(physical, true);
        ctx.stroke();
      }

      if (state.display.showElementIds) {
        const centroid = physical.reduce((acc, point) => ({
          x: acc.x + point.x / physical.length,
          y: acc.y + point.y / physical.length,
        }), { x: 0, y: 0 });
        const screen = worldToScreen(centroid);
        ctx.fillStyle = "rgba(224, 229, 246, 0.84)";
        ctx.font = "9px ui-monospace, monospace";
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(String(element.id), screen.x, screen.y);
      }
      elementIndex += 1;
    }
  }

  function pathPolygon(points, close) {
    ctx.beginPath();
    points.forEach((point, index) => {
      const screen = worldToScreen(point);
      if (index === 0) ctx.moveTo(screen.x, screen.y);
      else ctx.lineTo(screen.x, screen.y);
    });
    if (close || points.length >= 3) ctx.closePath();
  }

  function drawNodalContour(physical, values, nodeOrder) {
    if (physical.length < 3) return;
    const points = physical.map((point, index) => {
      const valueIndex = nodeOrder.has(point.id) ? nodeOrder.get(point.id) : index;
      return { x: point.x, y: point.y, value: values[valueIndex] || 0 };
    });

    for (let i = 1; i < points.length - 1; i += 1) {
      const triangle = [points[0], points[i], points[i + 1]];
      const centroid = triangle.reduce((acc, point) => ({
        x: acc.x + point.x / 3,
        y: acc.y + point.y / 3,
      }), { x: 0, y: 0 });
      const value = triangleValue(...triangle, centroid);
      ctx.fillStyle = normalizedColor(value, contourDomain());
      const screenTriangle = triangle.map(worldToScreen);
      ctx.beginPath();
      screenTriangle.forEach((point, index) => {
        if (index === 0) ctx.moveTo(point.x, point.y);
        else ctx.lineTo(point.x, point.y);
      });
      ctx.closePath();
      ctx.fill();
    }
  }

  function drawSelectedSets(nodes) {
    if (!state.model || !state.display.highlightSets || state.selectedSets.size === 0) return;
    ctx.save();
    for (const [name, kind] of state.selectedSets) {
      const set = kind === "nset" ? state.model.nsets.get(name) : state.model.elsets.get(name);
      if (!set) continue;
      ctx.fillStyle = "rgba(0, 211, 255, 0.18)";
      ctx.strokeStyle = "rgba(0, 211, 255, 0.88)";
      ctx.lineWidth = 1.7;
      if (kind === "nset") {
        for (const nodeId of set) {
          const node = nodes.get(nodeId);
          if (!node) continue;
          const screen = worldToScreen(node);
          ctx.beginPath();
          ctx.arc(screen.x, screen.y, 6, 0, Math.PI * 2);
          ctx.fill();
          ctx.stroke();
        }
      } else {
        for (const elementId of set) {
          const element = state.model.elements.get(elementId);
          if (!element) continue;
          const points = polygonFromConnectivity(element.connectivity, nodes);
          if (points.length < 2) continue;
          pathPolygon(points);
          if (points.length >= 3) ctx.fill();
          ctx.stroke();
        }
      }
    }
    ctx.restore();
  }

  function drawNodesAndIds(nodes) {
    if (!state.model) return;
    ctx.textAlign = "left";
    ctx.textBaseline = "alphabetic";
    for (const node of nodes.values()) {
      const screen = worldToScreen(node);
      if (state.display.showNodeDots) {
        ctx.fillStyle = "#f2f4fb";
        ctx.beginPath();
        ctx.arc(screen.x, screen.y, 2.1, 0, Math.PI * 2);
        ctx.fill();
      }
      if (state.display.showNodeIds) {
        ctx.fillStyle = "rgba(184, 203, 255, 0.92)";
        ctx.font = "9px ui-monospace, monospace";
        ctx.fillText(String(node.id), screen.x + 5, screen.y - 5);
      }
    }
  }

  function drawUndeformedOverlay() {
    if (!state.model || !state.results || !state.display.undeformedOverlay) return;
    if (deformationFactor() === 0) return;
    ctx.save();
    ctx.strokeStyle = "rgba(216, 198, 113, 0.62)";
    ctx.setLineDash([4, 3]);
    ctx.lineWidth = 1;
    for (const element of state.model.elements.values()) {
      const points = polygonFromConnectivity(element.connectivity, state.model.nodes);
      if (points.length < 2) continue;
      pathPolygon(points, true);
      ctx.stroke();
    }
    ctx.restore();
  }

  function drawBCsAndLoads(nodes) {
    if (!state.model) return;
    if (state.display.showBCs) {
      ctx.fillStyle = "#ff922b";
      ctx.strokeStyle = "rgba(255, 146, 43, 0.45)";
      for (const bc of state.model.bcs) {
        const ids = Array.isArray(bc && bc.nodes) ? bc.nodes : [];
        for (const nodeId of ids) {
          const node = nodes.get(Number(nodeId));
          if (!node) continue;
          const screen = worldToScreen(node);
          ctx.fillRect(screen.x - 4, screen.y - 4, 8, 8);
          ctx.strokeRect(screen.x - 4.5, screen.y - 4.5, 9, 9);
        }
      }
    }

    if (!state.display.showLoads) return;
    ctx.strokeStyle = "#ff4747";
    ctx.fillStyle = "#ff4747";
    ctx.lineWidth = 2;
    for (const load of state.model.cloads) {
      if (!Array.isArray(load) || load.length < 3) continue;
      const node = nodes.get(Number(load[0]));
      if (!node) continue;
      const dof = Number(load[1]);
      const magnitude = number(load[2]);
      const direction = dof === 1 ? Math.sign(magnitude) : 0;
      const vertical = dof === 2 ? -Math.sign(magnitude) : 0;
      const start = worldToScreen(node);
      const length = 26;
      const end = {
        x: start.x + direction * length,
        y: start.y + vertical * length,
      };
      if (!direction && !vertical) continue;
      ctx.beginPath();
      ctx.moveTo(start.x, start.y);
      ctx.lineTo(end.x, end.y);
      ctx.stroke();
      const angle = Math.atan2(end.y - start.y, end.x - start.x);
      ctx.beginPath();
      ctx.moveTo(end.x, end.y);
      ctx.lineTo(end.x - 7 * Math.cos(angle - 0.42), end.y - 7 * Math.sin(angle - 0.42));
      ctx.lineTo(end.x - 7 * Math.cos(angle + 0.42), end.y - 7 * Math.sin(angle + 0.42));
      ctx.closePath();
      ctx.fill();
    }
  }

  function drawLegend() {
    const domain = contourDomain();
    if (!domain) return;
    const { height } = canvasSize();
    const x = canvasSize().width - 57;
    const y = 22;
    const legendHeight = Math.min(220, Math.max(120, height - 100));
    const width = 17;

    ctx.save();
    ctx.fillStyle = "rgba(13, 14, 22, 0.82)";
    ctx.strokeStyle = "rgba(151, 162, 205, 0.5)";
    ctx.beginPath();
    ctx.roundRect(x - 6, y - 20, 62, legendHeight + 38, 3);
    ctx.fill();
    ctx.stroke();

    const steps = 60;
    for (let i = 0; i < steps; i += 1) {
      ctx.fillStyle = jetColor(i / (steps - 1));
      ctx.fillRect(x, y + legendHeight - (i + 1) * legendHeight / steps, width, legendHeight / steps + 0.5);
    }
    ctx.strokeStyle = "rgba(255,255,255,0.35)";
    ctx.strokeRect(x, y, width, legendHeight);

    ctx.fillStyle = "#dfe4f5";
    ctx.font = "10px ui-monospace, monospace";
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";
    ctx.fillText(domain.mode === "vm" ? "von Mises" : "|U|", x - 1, y - 11);
    ctx.fillText(domain.mode === "vm" ? "MPa" : "mm", x - 1, y + legendHeight + 9);

    const ticks = 6;
    for (let i = 0; i < ticks; i += 1) {
      const t = i / (ticks - 1);
      const tickY = y + legendHeight * (1 - t);
      const value = domain.min + (domain.max - domain.min) * t;
      const displayed = domain.mode === "vm" ? value / 1e6 : value * 1000;
      ctx.strokeStyle = "rgba(255,255,255,0.55)";
      ctx.beginPath();
      ctx.moveTo(x + width, tickY);
      ctx.lineTo(x + width + 3, tickY);
      ctx.stroke();
      ctx.fillStyle = "#e2e7f8";
      ctx.textAlign = "left";
      ctx.fillText(fmt(displayed, 3), x + width + 5, tickY);
    }
    ctx.restore();
  }

  function draw3DViewportDecorations() {
    const { width, height } = canvasSize();
    const originX = 14;
    const originY = height - 19;
    const axisLength = 37;
    const axes = [
      { axis: [1, 0, 0], color: "#ff7a7a", label: "X" },
      { axis: [0, 1, 0], color: "#76d37b", label: "Y" },
      { axis: [0, 0, 1], color: "#7aa2f4", label: "Z" },
    ];
    ctx.save();
    ctx.lineWidth = 1.5;
    ctx.font = "10px ui-monospace, monospace";
    ctx.textAlign = "left";
    ctx.textBaseline = "alphabetic";
    for (const item of axes) {
      const direction = matrixVectorMultiply(state.camera3D.rotation, item.axis);
      const end = {
        x: originX + direction[0] * axisLength,
        y: originY - direction[1] * axisLength,
      };
      ctx.strokeStyle = item.color;
      ctx.fillStyle = item.color;
      ctx.beginPath();
      ctx.moveTo(originX, originY);
      ctx.lineTo(end.x, end.y);
      ctx.stroke();
      ctx.fillText(item.label, end.x + 3, end.y + 3);
    }
    ctx.restore();

    const targetMeters = 90 / state.camera3D.scale;
    const exponent = Math.floor(Math.log10(Math.max(targetMeters, 1e-12)));
    const base = Math.pow(10, exponent);
    const multiplier = [1, 2, 5, 10].find((item) => item * base >= targetMeters) || 10;
    const barMeters = multiplier * base;
    const barPixels = barMeters * state.camera3D.scale;
    const barX = Math.max(58, width - barPixels - 22);
    const barY = height - 15;
    ctx.save();
    ctx.strokeStyle = "#e5e9f8";
    ctx.fillStyle = "#e5e9f8";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(barX, barY - 4);
    ctx.lineTo(barX, barY);
    ctx.lineTo(barX + barPixels, barY);
    ctx.lineTo(barX + barPixels, barY - 4);
    ctx.stroke();
    ctx.font = "10px ui-monospace, monospace";
    ctx.textAlign = "center";
    ctx.textBaseline = "alphabetic";
    ctx.fillText(fmtLength(barMeters), barX + barPixels / 2, barY - 7);
    ctx.restore();
  }

  function drawViewportDecorations() {
    const { width, height } = canvasSize();

    // Coordinate axes (model +X right, +Y up).
    const originX = 14;
    const originY = height - 19;
    const axisLength = 37;
    ctx.save();
    ctx.lineWidth = 1.5;
    ctx.strokeStyle = "#ff7a7a";
    ctx.beginPath();
    ctx.moveTo(originX, originY);
    ctx.lineTo(originX + axisLength, originY);
    ctx.stroke();
    ctx.strokeStyle = "#76d37b";
    ctx.beginPath();
    ctx.moveTo(originX, originY);
    ctx.lineTo(originX, originY - axisLength);
    ctx.stroke();
    ctx.fillStyle = "#d3d9ef";
    ctx.font = "10px ui-monospace, monospace";
    ctx.textAlign = "left";
    ctx.textBaseline = "alphabetic";
    ctx.fillText("X", originX + axisLength + 3, originY + 3);
    ctx.fillText("Y", originX - 2, originY - axisLength - 3);
    ctx.restore();

    // Scale bar: choose a round physical length near 90 px.
    const targetMeters = 90 / state.camera.scale;
    const exponent = Math.floor(Math.log10(Math.max(targetMeters, 1e-12)));
    const base = Math.pow(10, exponent);
    const multiplier = [1, 2, 5, 10].find((item) => item * base >= targetMeters) || 10;
    const barMeters = multiplier * base;
    const barPixels = barMeters * state.camera.scale;
    const barX = Math.max(58, width - barPixels - 22);
    const barY = height - 15;
    ctx.save();
    ctx.strokeStyle = "#e5e9f8";
    ctx.fillStyle = "#e5e9f8";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(barX, barY - 4);
    ctx.lineTo(barX, barY);
    ctx.lineTo(barX + barPixels, barY);
    ctx.lineTo(barX + barPixels, barY - 4);
    ctx.stroke();
    ctx.font = "10px ui-monospace, monospace";
    ctx.textAlign = "center";
    ctx.textBaseline = "alphabetic";
    ctx.fillText(fmtLength(barMeters), barX + barPixels / 2, barY - 7);
    ctx.restore();
  }

  function render() {
    const { width, height } = canvasSize();
    ctx.fillStyle = "#1a1a2e";
    ctx.fillRect(0, 0, width, height);

    if (!state.model) {
      ctx.fillStyle = "#969cb5";
      ctx.font = "13px Inter, system-ui, sans-serif";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText("No model loaded", width / 2, height / 2);
      drawViewportDecorations();
      return;
    }

    const deformed = Boolean(state.results) && state.display.deformationPercent > 0;
    const nodes = displayNodes();
    if (state.is3D) {
      draw3DUndeformedOverlay();
      const projected = draw3DMesh(nodes, deformed);
      draw3DSelectedSets(projected);
      draw3DBCsAndLoads(nodes, projected);
      draw3DNodesAndIds(projected);
      draw3DNodalContour(projected);
      drawLegend();
      draw3DViewportDecorations();
      return;
    }
    if (deformed) {
      drawMesh(nodes, true);
    } else {
      drawMesh(state.model.nodes, false);
    }
    drawUndeformedOverlay();
    drawSelectedSets(nodes);
    if (state.display.wireframe) {
      // Wireframe is drawn by drawMesh. Keep this branch explicit for module state.
    }
    drawBCsAndLoads(nodes);
    drawNodesAndIds(nodes);
    drawLegend();
    drawViewportDecorations();
  }

  /* UI event wiring */
  dom["fit-button"].addEventListener("click", () => {
    fitView();
    render();
  });

  dom["wireframe-button"].addEventListener("click", () => {
    state.display.wireframe = !state.display.wireframe;
    updateToolbarStates();
    render();
  });

  function setContour(mode) {
    state.display.contourMode = state.display.contourMode === mode ? "none" : mode;
    if (state.results && state.module !== "results") setModule("results", false);
    fitView();
    updateToolbarStates();
    updateInterface();
    render();
  }

  dom["contour-u-button"].addEventListener("click", () => setContour("u"));
  dom["contour-vm-button"].addEventListener("click", () => setContour("vm"));

  dom["undeformed-overlay-button"].addEventListener("click", () => {
    state.display.undeformedOverlay = !state.display.undeformedOverlay;
    syncDisplayInputs();
    updateToolbarStates();
    render();
  });

  dom["deform-slider"].addEventListener("input", () => {
    state.display.deformationPercent = Number(dom["deform-slider"].value);
    updateToolbarStates();
    renderResultChips();
    dom["viewport-mode-label"].textContent = viewportLabel();
    render();
  });

  const checkboxMap = [
    ["show-node-dots", "showNodeDots"],
    ["show-node-ids", "showNodeIds"],
    ["show-element-ids", "showElementIds"],
    ["highlight-sets", "highlightSets"],
    ["show-bcs", "showBCs"],
    ["show-loads", "showLoads"],
    ["show-undeformed-overlay", "undeformedOverlay"],
  ];
  checkboxMap.forEach(([id, key]) => {
    dom[id].addEventListener("change", () => {
      state.display[key] = dom[id].checked;
      render();
    });
  });

  function setModule(module, update = true) {
    state.module = module;
    if (module === "model") {
      state.display.contourMode = "none";
      state.display.deformationPercent = 0;
    } else if (module === "mesh") {
      state.display.contourMode = "none";
      state.display.showNodeDots = true;
      state.display.showNodeIds = true;
      state.display.showElementIds = true;
    } else if (module === "results") {
      if (state.display.contourMode === "none") state.display.contourMode = "vm";
      if (state.display.deformationPercent === 0) state.display.deformationPercent = 5;
    }
    syncDisplayInputs();
    if (update) {
      fitView();
      updateToolbarStates();
      updateInterface();
      render();
    }
  }

  function syncDisplayInputs() {
    dom["show-node-dots"].checked = state.display.showNodeDots;
    dom["show-node-ids"].checked = state.display.showNodeIds;
    dom["show-element-ids"].checked = state.display.showElementIds;
    dom["highlight-sets"].checked = state.display.highlightSets;
    dom["show-bcs"].checked = state.display.showBCs;
    dom["show-loads"].checked = state.display.showLoads;
    dom["show-undeformed-overlay"].checked = state.display.undeformedOverlay;
  }

  document.querySelectorAll(".module-button").forEach((button) => {
    button.addEventListener("click", () => setModule(button.dataset.module));
  });

  dom["job-form"].addEventListener("submit", async (event) => {
    event.preventDefault();
    await runDeck(dom["deck-path"].value);
  });

  dom["job-select"].addEventListener("change", () => {
    switchJob(dom["job-select"].value);
  });

  async function runDeck(path) {
    const deckPath = String(path || "").trim();
    if (!deckPath) {
      const payload = { ok: false, stage: "usage", error: "deck path is empty" };
      logJSON("error", payload, "POST run");
      return payload;
    }
    const body = { path: deckPath };
    logJSON("request", body, "POST /api/run");
    dom["run-button"].disabled = true;
    dom["run-button"].textContent = "Running";
    try {
      const payload = await requestJSON("/api/run", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
      state.runResponse = payload;
      logJSON(payload.ok === false ? "error" : "ok", payload, "POST run");
      if (payload.ok !== false && payload.job) {
        if (!state.jobs.includes(payload.job)) state.jobs.push(payload.job);
        ensureJobOption(payload.job);
        await Promise.all([
          loadModel(payload.job, { silent: true }),
          loadResults(payload.job, { silent: true }),
        ]);
        logJSON("ok", { model: state.model ? state.model.raw : null, results: state.results ? state.results.raw : null }, `loaded ${payload.job}`);
        setModule("results");
        setStatus(`Run complete: ${payload.job}`);
      } else {
        updateInterface();
      }
      return payload;
    } catch (error) {
      const payload = transportError(error);
      state.runResponse = payload;
      logJSON("error", payload, "POST run");
      setStatus("Run request failed");
      return payload;
    } finally {
      dom["run-button"].disabled = false;
      dom["run-button"].textContent = "Run";
    }
  }

  async function switchJob(job) {
    state.job = String(job);
    ensureJobOption(state.job);
    const modelLoaded = await loadModel(state.job, { silent: true });
    const resultsLoaded = await loadResults(state.job, { silent: true });
    const payload = {
      ok: modelLoaded,
      job: state.job,
      model: state.model ? state.model.raw : null,
      results: state.results ? state.results.raw : null,
      modelError: state.modelError,
      resultsError: state.resultsError,
    };
    logJSON(modelLoaded && resultsLoaded ? "ok" : "error", payload, `switch job ${state.job}`);
  }

  async function getAgent() {
    return requestJSON("/api/agent");
  }

  const commands = ["agent", "clear", "help", "model ", "results ", "run "];

  async function executeCommand(raw) {
    const command = String(raw || "").trim();
    if (!command) return null;
    state.commandHistory.push(command);
    state.historyIndex = state.commandHistory.length;
    logJSON("request", { command }, "console");
    const [verb, ...rest] = command.split(/\s+/);
    const argument = rest.join(" ");

    try {
      if (verb === "help") {
        const payload = {
          ok: true,
          commands: [
            { command: "run <path>", api: "POST /api/run" },
            { command: "model <job>", api: "GET /api/model?job=<job>" },
            { command: "results <job>", api: "GET /api/results/<job>" },
            { command: "agent", api: "GET /api/agent" },
            { command: "clear", api: null },
          ],
        };
        logJSON("info", payload, "help");
        return payload;
      }
      if (verb === "clear") {
        dom["console-output"].replaceChildren();
        state.consoleEntries = [];
        setStatus("Console cleared");
        return { ok: true, cleared: true };
      }
      if (verb === "run") {
        return await runDeck(argument);
      }
      if (verb === "model") {
        const job = argument || state.job;
        if (!await loadModel(job, { silent: true })) {
          logJSON("error", state.modelError, "console model");
          return state.modelError;
        }
        logJSON("ok", state.model.raw, `console model ${job}`);
        return state.model.raw;
      }
      if (verb === "results") {
        const job = argument || state.job;
        if (!await loadResults(job, { silent: true })) {
          logJSON("error", state.resultsError, "console results");
          return state.resultsError;
        }
        logJSON("ok", state.results.raw, `console results ${job}`);
        return state.results.raw;
      }
      if (verb === "agent") {
        const payload = await getAgent();
        logJSON("ok", payload, "GET agent");
        return payload;
      }
      const payload = {
        ok: false,
        stage: "usage",
        error: `unknown command: ${verb}`,
        hint: "try help",
      };
      logJSON("error", payload, "console");
      return payload;
    } catch (error) {
      const payload = transportError(error);
      logJSON("error", payload, "console");
      return payload;
    }
  }

  dom["console-form"].addEventListener("submit", async (event) => {
    event.preventDefault();
    const value = dom["command-input"].value;
    dom["command-input"].value = "";
    dom["command-hint"].textContent = "Tab completes a command · ↑/↓ recalls history";
    await executeCommand(value);
  });

  dom["command-input"].addEventListener("input", () => {
    const value = dom["command-input"].value.trimStart();
    const match = commands.find((command) => command.startsWith(value) && value && command !== value);
    dom["command-hint"].textContent = match
      ? `Tab → ${match}${match.endsWith(" ") ? "…" : ""}`
      : "Tab completes a command · ↑/↓ recalls history";
  });

  dom["command-input"].addEventListener("keydown", (event) => {
    if (event.key === "Tab") {
      event.preventDefault();
      const value = dom["command-input"].value.trimStart();
      const match = commands.find((command) => command.startsWith(value) && value && command !== value);
      if (match) {
        dom["command-input"].value = match;
        dom["command-hint"].textContent = match.endsWith(" ") ? "type a job or path…" : `Tab → ${match}`;
      }
    } else if (event.key === "ArrowUp") {
      event.preventDefault();
      if (!state.commandHistory.length) return;
      state.historyIndex = Math.max(0, state.historyIndex - 1);
      dom["command-input"].value = state.commandHistory[state.historyIndex] || "";
    } else if (event.key === "ArrowDown") {
      event.preventDefault();
      if (!state.commandHistory.length) return;
      state.historyIndex = Math.min(state.commandHistory.length, state.historyIndex + 1);
      dom["command-input"].value =
        state.commandHistory[state.historyIndex] || "";
    }
  });

  document.querySelectorAll("[data-command]").forEach((button) => {
    button.addEventListener("click", () => {
      dom["command-input"].value = button.dataset.command;
      dom["command-input"].focus();
    });
  });

  if (window.ResizeObserver) {
    new ResizeObserver(resizeCanvas).observe(canvas.parentElement);
  } else {
    window.addEventListener("resize", resizeCanvas);
  }

  if (window.__SOLVIX_EXPOSE_INTERNALS__) {
    window.SolvixUI3D = {
      HEX_EDGE_INDICES,
      TET_EDGE_INDICES,
      arcballVector,
      build3DEdgeRecords,
      edgeIndicesForElement,
      fitOrthographic,
      initial3DRotation,
      is3DModel,
      matrixIdentity,
      matrixMultiply,
      matrixTranspose,
      matrixVectorMultiply,
      painterSortEdges,
      panOrthographic,
      projectOrthographic,
      reconstructC3D8Depth,
      rotationBetweenUnitVectors,
      rotationXY,
      screenToCenterWorld,
      zoomOrthographicAt,
    };
  }

  async function initialize() {
    const params = new URLSearchParams(window.location.search);
    const requestedJob = params.get("job") || "cantilever3d";
    state.job = requestedJob;
    ensureJobOption(requestedJob);
    ensureJobOption("cantilever");
    resizeCanvas();
    logJSON("info", {
      ok: true,
      job: requestedJob,
      source: "startup",
      endpoints: ["/api/model", "/api/run", "/api/results/<job>", "/api/agent"],
      fallback: [
        "static/mock/model.json",
        "static/mock/results.json",
        "static/mock/model3d.json",
        "static/mock/results3d.json",
      ],
    }, "Solvix CAE");
    await loadModel(requestedJob, { silent: true });
    if (state.modelError) logJSON("info", state.modelError, "model source");
    if (state.model) logJSON("ok", state.model.raw, "GET model");
    await loadResults(requestedJob, { silent: true });
    if (state.resultsError) logJSON("info", state.resultsError, "results source");
    if (state.results) logJSON("ok", state.results.raw, "GET results");
    if (!state.results && requestedJob === "cantilever3d" && state.modelSource === "api") {
      logJSON("info", {
        ok: true,
        job: requestedJob,
        action: "auto-run",
        path: "examples/cantilever3d.inp",
      }, "startup");
      await runDeck("examples/cantilever3d.inp");
    }
    if (state.warnings.length) {
      logJSON("info", { ok: true, warnings: state.warnings }, "data consistency");
    }
    setModule(state.results ? "results" : "model");
  }

  initialize().catch((error) => {
    logJSON("error", transportError(error), "startup");
  });
}());
