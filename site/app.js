"use strict";

(function renderBenchmarks() {
  const benchmarks = Object.freeze([
    {
      name: "Cantilever family",
      release: "v0.1.0",
      value: 15.238,
      display: "15.238 mm",
      check: "Euler-Bernoulli tip-displacement target for CPS4 and C3D8 convergence",
    },
    {
      name: "Kirsch plate",
      release: "v0.1.0",
      value: 2.927,
      display: "2.927 MPa",
      check: "Extrapolated nodal von Mises stress against the 3.0 MPa analytical value",
    },
    {
      name: "First bending mode",
      release: "v0.3.0",
      value: 21.916,
      display: "21.916 Hz",
      check: "Cantilever modal frequency against the 20.886 Hz beam reference",
    },
    {
      name: "Critical load",
      release: "v0.3.0",
      value: 432081,
      display: "432081 N Pcr",
      check: "First linear buckling eigenvalue on the prestressed cantilever",
    },
    {
      name: "Thermal stress",
      release: "v0.4.0",
      value: 252.0,
      display: "252.0 MPa",
      check: "Sequential heat-to-stress benchmark through the coupled solve",
    },
    {
      name: "Dynamic amplification",
      release: "v0.4.0",
      value: 1.8,
      display: "1.8\u00d7 static",
      check: "Sudden-load transient peak relative to static equilibrium",
    },
    {
      name: "Nonlinear elastica",
      release: "v0.5.0",
      value: 0.9328,
      display: "0.9328",
      check: "Large-rotation cantilever ratio tracked toward the exact elastica",
    },
  ]);

  const tbody = document.querySelector("#benchmark-table tbody");
  const status = document.getElementById("benchmark-status");
  if (!tbody || !status) {
    return;
  }

  const rows = benchmarks.map((benchmark) => {
    const row = document.createElement("tr");
    const name = document.createElement("td");
    name.className = "benchmark-name";
    name.textContent = benchmark.name;

    const release = document.createElement("td");
    release.className = "benchmark-release";
    release.textContent = benchmark.release;

    const result = document.createElement("td");
    const resultValue = document.createElement("data");
    resultValue.className = "benchmark-value";
    resultValue.value = String(benchmark.value);
    resultValue.textContent = benchmark.display;
    result.appendChild(resultValue);

    const check = document.createElement("td");
    check.className = "benchmark-check";
    check.textContent = benchmark.check;

    row.append(name, release, result, check);
    return row;
  });

  tbody.replaceChildren(...rows);
  status.textContent = `Rendered ${benchmarks.length} static benchmark records.`;
  status.hidden = false;
})();
