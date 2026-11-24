'use strict';
/* =====================================================================
   Series RLC step response  |  L·i'' + R·i' + i/C = 0,  i(0)=0, vC(0)=0
   ===================================================================== */

/* ---------- 1. Mathematical core (pure functions, no DOM) ---------- */
function rlcModel(R, L, C, V) {
  const w0 = 1 / Math.sqrt(L * C);          // natural angular frequency  [rad/s]
  const alpha = R / (2 * L);                // damping factor (Neper freq.) [1/s]
  const zeta = alpha / w0;                  // damping ratio
  const Rc = 2 * Math.sqrt(L / C);          // critical resistance [Ohm]
  const m = { R, L, C, V, w0, alpha, zeta, Rc };
  m.type = R === 0 ? 'undamped' : Math.abs(zeta - 1) < 1e-6 ? 'critical' : zeta < 1 ? 'under' : 'over';

  if (m.type === 'under' || m.type === 'undamped') {
    const wd = Math.sqrt((w0 - alpha) * (w0 + alpha));
    m.wd = wd;
    m.i = t => V / (L * wd) * Math.exp(-alpha * t) * Math.sin(wd * t);
    m.v = t => V * (1 - Math.exp(-alpha * t) * (Math.cos(wd * t) + alpha / wd * Math.sin(wd * t)));
    m.tPkI = Math.atan2(wd, alpha) / wd;
    m.tPkV = Math.PI / wd;
    m.overshoot = Math.exp(-alpha * Math.PI / wd);
    m.vPk = V * (1 + m.overshoot);
  } else if (m.type === 'critical') {
    m.i = t => V / L * t * Math.exp(-alpha * t);
    m.v = t => V * (1 - (1 + alpha * t) * Math.exp(-alpha * t));
    m.tPkI = 1 / alpha;
    m.overshoot = 0; m.vPk = V; m.tPkV = null;
  } else {
    const beta = Math.sqrt((alpha - w0) * (alpha + w0));
    const s2 = -(alpha + beta), s1 = -(w0 * w0) / (alpha + beta);   // s1 = -alpha+beta, computed without cancellation
    m.beta = beta; m.s1 = s1; m.s2 = s2;
    m.i = t => V / (2 * beta * L) * Math.exp(s1 * t) * -Math.expm1(-2 * beta * t);   // (e^{s1 t}-e^{s2 t})
    m.v = t => V * (1 + (s2 * Math.exp(s1 * t) - s1 * Math.exp(s2 * t)) / (s1 - s2));
    m.tPkI = Math.log((alpha + beta) / w0) / beta;
    m.overshoot = 0; m.vPk = V; m.tPkV = null;
  }
  m.iPk = m.i(m.tPkI);

  // time window long enough for the response to settle
  if (m.type === 'undamped') m.tEnd = 6 * 2 * Math.PI / w0;
  else if (m.type === 'over') m.tEnd = 7 / -m.s1;
  else if (m.type === 'under') {
    const per = 2 * Math.PI / m.wd;
    m.tEnd = Math.max(7 / alpha, m.overshoot > 1e-3 ? 1.1 * m.tPkV : 0);   // only stretch to the overshoot peak if it is visible
    if (m.tEnd > 40 * per) { m.tEnd = Math.max(40 * per, 1.1 * m.tPkV); m.truncated = true; }   // keep oscillations resolvable
  } else m.tEnd = 7 / alpha;

  // 2 % settling time of vC, found from the exact solution
  m.settle = Infinity;
  if (m.type === 'under') {
    // |vC-V| <= V e^{-at} sqrt(1+(a/wd)^2): in band after tEnv; scan back from there to the last excursion
    const per = 2 * Math.PI / m.wd, dt = per / 200;
    let t = Math.log(50 * Math.sqrt(1 + (alpha / m.wd) ** 2)) / alpha;
    while (t > 0 && Math.abs(m.v(t) - V) <= 0.02 * V) t -= dt;
    m.settle = Math.max(t + dt, 0);
  } else if (m.type !== 'undamped') {
    const N = 20000; let last = -1;
    for (let k = 0; k <= N; k++) if (Math.abs(m.v(m.tEnd * k / N) - V) > 0.02 * V) last = k;
    m.settle = last >= N ? Infinity : m.tEnd * (last + 1) / N;
  }
  return m;
}

// sample the exact solution; uniform spacing for oscillatory cases, denser near t=0 otherwise
function sampleModel(m, tEnd = m.tEnd, forceUniform = false) {
  const osc = m.type === 'under' || m.type === 'undamped';
  const N = osc || forceUniform ? Math.min(6000, Math.max(1000, Math.ceil(tEnd / (2 * Math.PI / (m.wd || m.w0)) * 50))) : 1500;
  const t = [], i = [], v = [];
  for (let k = 0; k <= N; k++) {
    const x = k / N, tk = osc || forceUniform ? tEnd * x : tEnd * x * x;
    t.push(tk); i.push(m.i(tk)); v.push(m.v(tk));
  }
  return { t, i, v };
}

if (typeof module !== 'undefined' && module.exports) module.exports = { rlcModel, sampleModel };

/* ---------- 2. Formatting helpers ---------- */
const PFX = [[1e9, 'G'], [1e6, 'M'], [1e3, 'k'], [1, ''], [1e-3, 'm'], [1e-6, 'µ'], [1e-9, 'n'], [1e-12, 'p'], [1e-15, 'f']];
function si(x, u, d = 4) {
  if (!isFinite(x)) return (x > 0 ? '∞ ' : '−∞ ') + u;
  if (x === 0) return '0 ' + u;
  const a = Math.abs(x), p = PFX.find(q => a >= q[0] * 0.99995) || PFX[PFX.length - 1];
  return String(+(x / p[0]).toPrecision(d)) + ' ' + p[1] + u;
}
function n(x, d = 4) {            // number as HTML, scientific for very large/small
  if (!isFinite(x)) return x > 0 ? '∞' : '−∞';
  if (x === 0) return '0';
  const a = Math.abs(x);
  if (a >= 1e5 || a < 1e-3) { const [mm, e] = x.toExponential(d - 1).split('e'); return `${+mm}×10<sup>${+e}</sup>`; }
  return String(+x.toPrecision(d));
}
const res = s => `<span class="res">${s}</span>`;
function timeUnit(t) { return t < 1e-6 ? { n: 'ns', f: 1e9 } : t < 1e-3 ? { n: 'µs', f: 1e6 } : t < 1 ? { n: 'ms', f: 1e3 } : { n: 's', f: 1 }; }
function pickUnit(max, b) { return max >= 1 ? { n: b, f: 1 } : max >= 1e-3 ? { n: 'm' + b, f: 1e3 } : max >= 1e-6 ? { n: 'µ' + b, f: 1e6 } : { n: 'n' + b, f: 1e9 }; }
const TYPE_NAME = { under: 'Underdamped', critical: 'Critically damped', over: 'Overdamped', undamped: 'Undamped (R = 0)' };
const TYPE_CLASS = { under: 'b-under', critical: 'b-crit', over: 'b-over', undamped: 'b-und' };
const COLOR = { under: '#27d3ff', critical: '#3ddc97', over: '#ffb84a', user: '#e6eefc' };

/* ---------- 3. Input validation ---------- */
const LIM = {
  R: { name: 'Resistance', unit: 'Ω', min: 0, max: 1e7, zero: true },
  L: { name: 'Inductance', unit: 'H', min: 1e-9, max: 1e4 },
  C: { name: 'Capacitance', unit: 'µF', min: 1e-3, max: 1e7 },
  V: { name: 'Input voltage', unit: 'V', min: 1e-6, max: 1e6 }
};
function parseField(k, raw) {
  const q = LIM[k], s = String(raw).trim();
  if (s === '') return { err: `${q.name} is required.` };
  if (!/^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/.test(s)) return { err: 'Enter a plain number such as 0.01 or 1e-3 (no letters or units).' };
  const x = Number(s);
  if (!isFinite(x)) return { err: 'That number is too large to use.' };
  if (x < 0) return { err: `${q.name} cannot be negative.` };
  if (x === 0 && !q.zero) return { err: `${q.name} must be greater than 0 (a zero value makes the circuit equation undefined).` };
  if (x !== 0 && x < q.min) return { err: `Too small. Minimum supported is ${q.min} ${q.unit}.` };
  if (x > q.max) return { err: `Too large. Maximum supported is ${q.max} ${q.unit}.` };
  return { val: x };
}

/* ---------- 4. User interface ---------- */
const PRESETS = [
  { name: 'Underdamped', R: 5, L: 0.01, C: 100, V: 12 },
  { name: 'Critically damped', R: 20, L: 0.01, C: 100, V: 12 },
  { name: 'Overdamped', R: 100, L: 0.01, C: 100, V: 12 },
  { name: 'PSU output filter', R: 0.2, L: 0.000047, C: 470, V: 5 },
  { name: 'Lossless (R = 0)', R: 0, L: 0.01, C: 100, V: 12 }
];
const DEFAULT = PRESETS[0];

function initUI() {
  const $ = id => document.getElementById(id);
  const MUTED = '#8fa3c7', GRID = 'rgba(110,160,255,.12)';
  Chart.defaults.color = MUTED;
  Chart.defaults.font.family = "'Sora','Segoe UI',system-ui,sans-serif";
  let cur = null;

  /* charts */
  const opts = legend => ({
    responsive: true, maintainAspectRatio: false, parsing: false, normalized: true, animation: { duration: 650 },
    interaction: { mode: 'nearest', axis: 'x', intersect: false },
    plugins: {
      legend: { display: legend, labels: { color: '#cfe0ff', usePointStyle: true, boxWidth: 8 } },
      tooltip: {
        backgroundColor: 'rgba(8,16,34,.96)', borderColor: 'rgba(110,160,255,.3)', borderWidth: 1,
        callbacks: {
          title: it => 't = ' + (+it[0].parsed.x.toPrecision(4)) + ' ' + it[0].chart.$u.t,
          label: c => c.dataset.label + ': ' + (+c.parsed.y.toPrecision(4)) + ' ' + c.chart.$u.y
        }
      }
    },
    scales: {
      x: { type: 'linear', min: 0, title: { display: true, text: 'Time', color: MUTED }, grid: { color: GRID }, ticks: { maxTicksLimit: 9 } },
      y: { title: { display: true, text: '', color: MUTED }, grid: { color: GRID } }
    }
  });
  const mk = (id, legend) => new Chart($(id), { type: 'line', data: { datasets: [] }, options: opts(legend) });
  const chI = mk('chartI', false), chV = mk('chartV', false), cI = mk('cmpI', true), cV = mk('cmpV', true), hero = mk('heroChart', false);
  hero.options.plugins.tooltip.enabled = false;

  const pts = (S, key, tf, yf) => S.t.map((t, k) => ({ x: t * tf, y: S[key][k] * yf }));
  const line = (label, data, color, x = {}) => Object.assign({ label, data, borderColor: color, backgroundColor: color + '14', borderWidth: 2.4, pointRadius: 0, pointHoverRadius: 4, tension: 0 }, x);
  const dot = (label, data, color, r = 6) => ({ label, type: 'scatter', data, backgroundColor: color, borderColor: '#fff', borderWidth: 1.5, pointRadius: r, pointHoverRadius: r + 2 });
  function setChart(ch, datasets, tEnd, tu, yu, ylab, anim) {
    ch.$u = { t: tu.n, y: yu.n }; ch.data.datasets = datasets;
    ch.options.scales.x.max = tEnd * tu.f; ch.options.scales.x.title.text = `Time (${tu.n})`;
    ch.options.scales.y.title.text = `${ylab} (${yu.n})`; ch.update(anim ? undefined : 'none');
  }

  /* hero preview */
  (function () {
    const m = rlcModel(5, 0.01, 100e-6, 12), S = sampleModel(m), tu = timeUnit(m.tEnd);
    setChart(hero, [line('v_C(t)', pts(S, 'v', tu.f, 1), COLOR.under, { fill: true }),
      line('V', [{ x: 0, y: 12 }, { x: m.tEnd * tu.f, y: 12 }], '#3ddc97', { borderDash: [6, 6], borderWidth: 1.2, backgroundColor: 'transparent' })], m.tEnd, tu, { n: 'V' }, 'Voltage', true);
  })();

  /* validation UI */
  const inputs = { R: $('R'), L: $('L'), C: $('C'), V: $('V') };
  function showAlert(msg, note) { const a = $('alert'); a.textContent = msg; a.className = 'alert show' + (note ? ' note' : ''); }
  function hideAlert() { $('alert').className = 'alert'; }
  function validate() {
    const vals = {}; let ok = true, first = null;
    for (const k in inputs) {
      const r = parseField(k, inputs[k].value), box = inputs[k].closest('.field');
      box.classList.toggle('err', !!r.err); box.querySelector('.msg').textContent = r.err || '';
      if (r.err) { ok = false; first = first || k; } else vals[k] = r.val;
    }
    if (!ok) showAlert('Fix the highlighted input(s). Results below still show the last valid simulation.');
    return ok ? vals : null;
  }

  /* run everything */
  function run(anim) {
    const vals = validate(); if (!vals) return false;
    const m = rlcModel(vals.R, vals.L, vals.C * 1e-6, vals.V); cur = m;
    hideAlert();
    if (m.type === 'undamped') showAlert('R = 0 is the ideal lossless case: nothing dissipates energy, so the response oscillates forever and never settles.', true);
    else if (m.zeta > 1e3) showAlert('Extreme damping (ζ > 1000): the response is extremely slow compared with its initial rise. Hover the graph to read early values.', true);
    else if (m.truncated) showAlert('Very light damping: the response rings for hundreds of cycles. The graphs show the first 40 cycles; the settling time shown is still exact.', true);
    renderDiagram(m); renderResults(m); renderSteps(m); renderGraphs(m, anim); renderCompare(m, anim);
    return true;
  }

  function renderDiagram(m) {
    $('lblR').textContent = 'R = ' + si(m.R, 'Ω'); $('lblL').textContent = 'L = ' + si(m.L, 'H');
    $('lblC').textContent = 'C = ' + si(m.C, 'F'); $('lblV').textContent = 'V = ' + si(m.V, 'V');
  }

  function renderResults(m) {
    const stat = (l, v) => `<div class="glass stat"><span>${l}</span><strong>${v}</strong></div>`;
    const second = m.type === 'over' ? stat('Real-root spread β', si(m.beta, 's⁻¹')) : m.type === 'critical' ? stat('Repeated root s', si(-m.alpha, 's⁻¹')) : stat('Damped frequency ω<sub>d</sub>', si(m.wd, 'rad/s'));
    $('results').innerHTML =
      stat('Natural frequency ω<sub>n</sub>', si(m.w0, 'rad/s')) + stat('Damping factor α', si(m.alpha, 's⁻¹')) +
      stat('Damping ratio ζ', n(m.zeta)) + stat('Critical resistance R<sub>c</sub>', si(m.Rc, 'Ω')) +
      `<div class="glass stat"><span>Classification</span><strong><span class="badge ${TYPE_CLASS[m.type]}">${TYPE_NAME[m.type]}</span></strong></div>` + second +
      stat('Peak current', `${si(m.iPk, 'A')}<br><small>at ${si(m.tPkI, 's')}</small>`) +
      stat('Peak capacitor voltage', m.type === 'under' || m.type === 'undamped' ? `${si(m.vPk, 'V')}<br><small>at ${si(m.tPkV, 's')}</small>` : `${si(m.V, 'V')}<br><small>final value, no overshoot</small>`) +
      stat('Overshoot', (m.overshoot * 100 < 0.005 ? '0' : (+(m.overshoot * 100).toPrecision(4))) + ' %') +
      stat('Settling time (±2 %)', isFinite(m.settle) ? si(m.settle, 's') : 'never (no damping)') +
      stat('Final capacitor voltage', si(m.V, 'V')) + stat('Final current', '0 A');
    const f = m.type === 'under' || m.type === 'undamped' ? m.wd / (2 * Math.PI) : 0;
    const head = `<h4>What this means</h4>`;
    const txt = {
      under: `Because ζ = ${n(m.zeta)} is below 1, the resistor is too small to stop the inductor and capacitor from trading energy back and forth. The current rises from zero to ${si(m.iPk, 'A')} at ${si(m.tPkI, 's')}, then swings back and forth about ${si(f, 'Hz')} while fading. The capacitor voltage shoots past the supply, peaking at ${si(m.vPk, 'V')} (${(+(m.overshoot * 100).toPrecision(3))} % overshoot), then settles to ${si(m.V, 'V')} after about ${si(m.settle, 's')}. In a power supply, that overshoot is what can stress downstream parts.`,
      critical: `R equals the critical resistance, so the circuit returns to rest as fast as possible without oscillating. The current climbs to ${si(m.iPk, 'A')} at ${si(m.tPkI, 's')} and then dies away smoothly. The capacitor charges to ${si(m.V, 'V')} with no overshoot and stays within 2 % after ${si(m.settle, 's')}.`,
      over: `Because ζ = ${n(m.zeta)} is above 1, the resistor dominates. The current is held down to a peak of ${si(m.iPk, 'A')} at ${si(m.tPkI, 's')} and then falls slowly. The capacitor charges without any overshoot, but it creeps toward ${si(m.V, 'V')} and needs about ${si(m.settle, 's')} to get within 2 %. The circuit is safe but sluggish.`,
      undamped: `With R = 0 nothing removes energy, so the circuit is an ideal LC oscillator. The current swings between 0 and ±${si(m.iPk, 'A')} forever, and the capacitor voltage swings between 0 and exactly 2V = ${si(2 * m.V, 'V')}. Real circuits always have some resistance, which is why real responses decay.`
    }[m.type];
    $('meaning').innerHTML = head + `<p>${txt}</p>`;
  }

  function renderSteps(m) {
    const { R, L, C, V, w0, alpha, zeta, Rc } = m;
    const st = [];
    const add = (t, f, p) => st.push(`<li><b>${t}</b>${p ? p : ''}${f ? `<div class="formula">${f}</div>` : ''}</li>`);
    add('Convert units to SI', `C = ${n(C * 1e6)} µF × 10<sup>−6</sup> = ${res(n(C) + ' F')}`, `R = ${n(R)} Ω, L = ${n(L)} H and V = ${n(V)} V are already in SI units.`);
    add('Natural angular frequency', `ω<sub>n</sub> = 1/√(L·C) = 1/√(${n(L)} × ${n(C)}) = ${res(n(w0) + ' rad/s')}`, `That is ${si(w0 / (2 * Math.PI), 'Hz')}.`);
    add('Damping factor', `α = R/(2L) = ${n(R)}/(2 × ${n(L)}) = ${res(n(alpha) + ' s⁻¹')}`);
    add('Damping ratio', `ζ = α/ω<sub>n</sub> = ${n(alpha)}/${n(w0)} = ${res(n(zeta))}`);
    add('Critical resistance', `R<sub>c</sub> = 2·√(L/C) = 2·√(${n(L)}/${n(C)}) = ${res(n(Rc) + ' Ω')}`, `R is ${R < Rc ? 'below' : R > Rc ? 'above' : 'equal to'} R<sub>c</sub> (R/R<sub>c</sub> = ${n(R / Rc)}), which matches ζ.`);
    add('Classification', `${m.type === 'under' || m.type === 'undamped' ? 'ζ < 1' : m.type === 'over' ? 'ζ > 1' : 'ζ = 1'} → ${res(TYPE_NAME[m.type].toLowerCase())}`);
    if (m.type === 'under' || m.type === 'undamped') {
      add('Damped frequency and solution',
        `ω<sub>d</sub> = √(ω<sub>n</sub>² − α²) = √(${n(w0)}² − ${n(alpha)}²) = ${res(n(m.wd) + ' rad/s')}<br>i(t) = ${n(V / (L * m.wd))} · e<sup>−${n(alpha)}t</sup> · sin(${n(m.wd)}t) &nbsp;A<br>v<sub>C</sub>(t) = ${n(V)} · [1 − e<sup>−${n(alpha)}t</sup>(cos(${n(m.wd)}t) + ${n(alpha / m.wd)}·sin(${n(m.wd)}t))] &nbsp;V`);
    } else if (m.type === 'critical') {
      add('Solution (repeated root s = −α)', `i(t) = (V/L)·t·e<sup>−αt</sup> = ${n(V / L)} · t · e<sup>−${n(alpha)}t</sup> &nbsp;A<br>v<sub>C</sub>(t) = ${n(V)} · [1 − (1 + ${n(alpha)}t)·e<sup>−${n(alpha)}t</sup>] &nbsp;V`);
    } else {
      add('Roots and solution', `β = √(α² − ω<sub>n</sub>²) = √(${n(alpha)}² − ${n(w0)}²) = ${res(n(m.beta) + ' s⁻¹')}<br>s₁ = −α + β = ${res(n(m.s1) + ' s⁻¹')}, &nbsp; s₂ = −α − β = ${res(n(m.s2) + ' s⁻¹')}<br>i(t) = ${n(V / (2 * m.beta * L))} · (e<sup>${n(m.s1)}t</sup> − e<sup>${n(m.s2)}t</sup>) &nbsp;A<br>v<sub>C</sub>(t) = ${n(V)} · [1 + (s₂e<sup>s₁t</sup> − s₁e<sup>s₂t</sup>)/(s₁ − s₂)] &nbsp;V`);
    }
    add('Check the initial conditions', `i(0) = ${n(m.i(0))} A &nbsp; v<sub>C</sub>(0) = ${n(m.v(0))} V &nbsp; di/dt(0⁺) = V/L = ${n(V / L)} A/s`, 'The inductor current starts at zero and the capacitor starts uncharged, so the whole source voltage first appears across the inductor.');
    if (m.type === 'under' || m.type === 'undamped')
      add('Peak current', `t<sub>pk</sub> = atan(ω<sub>d</sub>/α)/ω<sub>d</sub> = ${res(si(m.tPkI, 's'))}<br>i<sub>max</sub> = i(t<sub>pk</sub>) = ${res(si(m.iPk, 'A'))}`);
    else if (m.type === 'critical')
      add('Peak current', `t<sub>pk</sub> = 1/α = ${res(si(m.tPkI, 's'))}<br>i<sub>max</sub> = 2V/(e·R) = ${res(si(m.iPk, 'A'))}`);
    else add('Peak current', `t<sub>pk</sub> = ln((α + β)/ω<sub>n</sub>)/β = ${res(si(m.tPkI, 's'))}<br>i<sub>max</sub> = i(t<sub>pk</sub>) = ${res(si(m.iPk, 'A'))}`);
    if (m.type === 'under' || m.type === 'undamped')
      add('Peak capacitor voltage', `t = π/ω<sub>d</sub> = ${res(si(m.tPkV, 's'))}<br>v<sub>C,max</sub> = V·(1 + e<sup>−απ/ω<sub>d</sub></sup>) = ${n(V)} × (1 + ${n(m.overshoot)}) = ${res(si(m.vPk, 'V'))}`, `Overshoot = ${res((+(m.overshoot * 100).toPrecision(4)) + ' %')}.`);
    else add('Peak capacitor voltage', `v<sub>C</sub> rises monotonically to V, so v<sub>C,max</sub> = ${res(si(V, 'V'))} (no overshoot)`);
    add('Settling time (±2 % of V)', isFinite(m.settle) ? `T<sub>s</sub> = ${res(si(m.settle, 's'))}${m.type === 'under' ? `&nbsp; (envelope estimate 4/α = ${si(4 / alpha, 's')})` : ''}` : `${res('does not settle')} (R = 0)`, 'Found as the last time |v<sub>C</sub> − V| exceeds 0.02·V, using the exact solution.');
    $('steps').innerHTML = st.join('');
  }

  /* single-simulation graphs */
  function renderGraphs(m, anim) {
    const S = sampleModel(m), tu = timeUnit(m.tEnd), iu = pickUnit(m.iPk, 'A'), vu = pickUnit(Math.max(m.vPk, m.V), 'V');
    const col = COLOR[m.type];
    $('legI').innerHTML = `<span><i style="background:${col}"></i>i(t)</span><span><i style="background:#ffb84a"></i>peak</span>`;
    $('legV').innerHTML = `<span><i style="background:${col}"></i>v<sub>C</sub>(t)</span><span><i style="background:#3ddc97"></i>final V, settling</span>`;
    const iPk = dot('Peak current', [{ x: m.tPkI * tu.f, y: m.iPk * iu.f }], '#ffb84a');
    setChart(chI, [line('i(t)', pts(S, 'i', tu.f, iu.f), col, { fill: true }), iPk, dot('Cursor', [{ x: 0, y: 0 }], '#fff', 5)], m.tEnd, tu, iu, 'Current', anim);
    const vds = [line('v_C(t)', pts(S, 'v', tu.f, vu.f), col, { fill: true }),
      line('Final value V', [{ x: 0, y: m.V * vu.f }, { x: m.tEnd * tu.f, y: m.V * vu.f }], '#3ddc97', { borderDash: [6, 6], borderWidth: 1.3, backgroundColor: 'transparent' }),
      dot('Peak v_C', m.tPkV ? [{ x: m.tPkV * tu.f, y: m.vPk * vu.f }] : [], '#ffb84a'),
      dot('Settling (2 %)', isFinite(m.settle) ? [{ x: m.settle * tu.f, y: m.v(m.settle) * vu.f }] : [], '#3ddc97', 5),
      dot('Cursor', [{ x: 0, y: 0 }], '#fff', 5)];
    setChart(chV, vds, m.tEnd, tu, vu, 'Voltage', anim);
    chI.$f = { tu, iu }; chV.$f = { tu, vu };
    updateCursor();
  }
  function updateCursor() {
    if (!cur) return;
    const m = cur, t = m.tEnd * $('slider').value / 1000, i = m.i(t), v = m.v(t), tu = chI.$f.tu, iu = chI.$f.iu, vu = chV.$f.vu;
    chI.data.datasets[2].data = [{ x: t * tu.f, y: i * iu.f }]; chV.data.datasets[4].data = [{ x: t * tu.f, y: v * vu.f }];
    chI.update('none'); chV.update('none');
    $('cursorRead').innerHTML = `t = <b>${si(t, 's')}</b> &nbsp; i = <b>${si(i, 'A')}</b> &nbsp; v<sub>C</sub> = <b>${si(v, 'V')}</b> &nbsp; v<sub>L</sub> + v<sub>R</sub> = <b>${si(m.V - v, 'V')}</b>`;
  }

  /* comparison: three separate simulations at ζ = 0.2, 1, 2.5 plus the user's R */
  function renderCompare(m, anim) {
    const Rc = m.Rc;
    const cases = [['under', 0.2], ['critical', 1], ['over', 2.5]].map(([k, z]) => ({ key: k, label: `${TYPE_NAME[k]} (ζ = ${z})`, m: rlcModel(z * Rc, m.L, m.C, m.V) }));
    cases.push({ key: 'user', label: `Your R = ${si(m.R, 'Ω')} (ζ = ${+m.zeta.toPrecision(3)})`, m });
    const tEnd = Math.max(...cases.slice(0, 3).map(c => c.m.tEnd));
    const tu = timeUnit(tEnd), iu = pickUnit(Math.max(...cases.map(c => c.m.iPk)), 'A'), vu = pickUnit(Math.max(...cases.map(c => c.m.vPk)), 'V');
    const mkds = (kv, f) => cases.map(c => {
      const S = sampleModel(c.m, tEnd, true);
      return line(c.label, pts(S, kv, tu.f, f), COLOR[c.key], c.key === 'user' ? { borderDash: [7, 5], borderWidth: 2, backgroundColor: 'transparent' } : { backgroundColor: 'transparent' });
    });
    setChart(cI, mkds('i', iu.f), tEnd, tu, iu, 'Current', anim);
    setChart(cV, mkds('v', vu.f), tEnd, tu, vu, 'Voltage', anim);
    const row = c => `<tr><td style="font-family:inherit;color:${COLOR[c.key]}">${c.label}</td><td>${si(c.m.R, 'Ω')}</td><td>${+c.m.zeta.toPrecision(3)}</td><td>${si(c.m.iPk, 'A')} @ ${si(c.m.tPkI, 's')}</td><td>${si(c.m.vPk, 'V')}</td><td>${+(c.m.overshoot * 100).toPrecision(3)} %</td><td>${isFinite(c.m.settle) ? si(c.m.settle, 's') : 'never'}</td></tr>`;
    $('cmpTable').innerHTML = '<tr><th>Case</th><th>R</th><th>ζ</th><th>Peak current</th><th>Peak v<sub>C</sub></th><th>Overshoot</th><th>Settling (±2 %)</th></tr>' + cases.map(row).join('');
  }

  /* wiring */
  $('presets').innerHTML = PRESETS.map((p, k) => `<button type="button" class="btn ghost sm" data-p="${k}">${p.name}</button>`).join('');
  function fill(p) { inputs.R.value = p.R; inputs.L.value = p.L; inputs.C.value = p.C; inputs.V.value = p.V; $('slider').value = 0; run(true); }
  $('presets').addEventListener('click', e => { const b = e.target.closest('[data-p]'); if (b) fill(PRESETS[b.dataset.p]); });
  $('form').addEventListener('submit', e => { e.preventDefault(); $('slider').value = 0; run(true); });
  $('resetBtn').addEventListener('click', () => fill(DEFAULT));
  $('cmpBtn').addEventListener('click', () => { if (run(true)) $('compare').scrollIntoView({ behavior: 'smooth' }); });
  $('cmpRun').addEventListener('click', () => run(true));
  $('slider').addEventListener('input', updateCursor);
  let timer; Object.values(inputs).forEach(el => el.addEventListener('input', () => { clearTimeout(timer); timer = setTimeout(() => run(false), 220); }));
  $('startBtn').addEventListener('click', () => setTimeout(() => inputs.R.focus({ preventScroll: true }), 700));

  /* navigation */
  const links = $('links');
  $('burger').addEventListener('click', () => { const o = links.classList.toggle('open'); $('burger').setAttribute('aria-expanded', o); });
  links.addEventListener('click', e => { if (e.target.tagName === 'A') { links.classList.remove('open'); $('burger').setAttribute('aria-expanded', false); } });
  const io = new IntersectionObserver(es => es.forEach(en => {
    if (en.isIntersecting) links.querySelectorAll('a').forEach(a => a.classList.toggle('active', a.getAttribute('href') === '#' + en.target.id));
  }), { rootMargin: '-40% 0px -55% 0px' });
  document.querySelectorAll('main section').forEach(s => io.observe(s));

  run(true);
}
if (typeof window !== 'undefined' && typeof document !== 'undefined') initUI();
