const $ = i => document.getElementById(i);

// ===== Data soal =====
const PR = [
  {n:'Soal 1', dim:1, f:'r - x^2', p:{r:[1,-3,4]}, d:[-4,4],
   t:'<b>Soal 1: x′ = r − x²</b>Buku 1 (Strogatz), hlm. 48, Contoh 3.1.1. Geser r: dua titik tetap untuk r > 0, satu untuk r = 0, tidak ada untuk r < 0.'},
  {n:'Soal 2', dim:2, f:'y', g:'-(g/l)*sin(x) - (k/m)*y', p:{g:[9.8,1,20], l:[1,.2,3], k:[.5,0,3], m:[1,.2,3]}, d:[-4,4,-4,4],
   t:'<b>Soal 2: Persamaan pendulum</b>Buku 3 (Khalil), hlm. 5–6. x₁ = θ, x₂ = θ′. Titik (0,0) menggantung ke bawah, titik (±π,0) terbalik ke atas.'},
  {n:'Soal 3', dim:2, f:'sin(pi*y) + (x-1)^2', g:'y^2 - y', p:{}, d:[-1,3,-1,2],
   t:'<b>Soal 3: x′ = sin(πy) + (x−1)², y′ = y² − y</b>Buku 7 (Notes on Diffy Qs), hlm. 355–356. Titik kritis (1,0) dan (1,1).'},
  {n:'Soal 4', dim:2, f:'x*(3 - x - 2*y)', g:'y*(2 - x - y)', p:{}, d:[-.5,4,-.5,3],
   t:'<b>Soal 4: Dua spesies berkompetisi</b>Buku 4, hlm. 170. x′ = x(3 − x − 2y), y′ = y(2 − x − y).'},
  {n:'Soal 5', dim:2, f:'x*(-6 - 2*x - 3*y)', g:'y*(1 - x - y)', p:{}, d:[-5,11,-10,3],
   t:'<b>Soal 5: Linearisasi sistem non-linier</b>Buku 6 (Judson), Bab 5.1, hlm. 288–300. x′ = x(−6 − 2x − 3y), y′ = y(1 − x − y).'},
  {n:'Sistem sendiri', dim:2, f:'x - x^3', g:'-y', p:{}, d:[-3,3,-3,3],
   t:'<b>Sistem sendiri</b>Ubah persamaan f dan g. Contoh awal: x′ = x − x³, y′ = −y.'}
];

let cur = 0, F, G, eqs = [], dom = [0,1,0,1], dim = 2;

// ===== Fungsi bantu =====
// Format angka; kelipatan π ditampilkan sebagai π
function fmt(v){
  if (Math.abs(v) < 5e-7) return '0';
  const q = v / Math.PI;
  if (Math.abs(q - Math.round(q)) < 1e-5){
    const r = Math.round(q);
    return (r === 1 ? '' : r === -1 ? '-' : r) + 'π';
  }
  return String(+v.toFixed(4));
}
const r6 = v => { v = Math.round(v * 1e6) / 1e6; return v === 0 ? 0 : v; };

// Mengubah teks persamaan menjadi fungsi JavaScript
function comp(e, pn, pv){
  e = e.replace(/\^/g, '**').replace(/\bpi\b/gi, 'PI');
  const f = new Function(...pn, 'x', 'y', 'with(Math){return (' + e + ')}');
  return (x, y) => f(...pv, x, y);
}

// ===== Jacobian (beda hingga terpusat) =====
function jac(a, b, x, y){
  const h = 1e-6;
  return [
    [(a(x+h,y) - a(x-h,y)) / (2*h), (a(x,y+h) - a(x,y-h)) / (2*h)],
    [(b(x+h,y) - b(x-h,y)) / (2*h), (b(x,y+h) - b(x,y-h)) / (2*h)]
  ];
}

// ===== Newton-Raphson untuk sistem 2 persamaan =====
function newton(a, b, x, y){
  for (let i = 0; i < 150; i++){
    const u = a(x,y), v = b(x,y);
    if (Math.hypot(u,v) < 1e-11) return [x,y];
    const J = jac(a,b,x,y), d = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    if (Math.abs(d) < 1e-14) return null;
    x -= (u*J[1][1] - v*J[0][1]) / d;
    y -= (J[0][0]*v - J[1][0]*u) / d;
    if (!isFinite(x+y) || Math.abs(x) > 1e4 || Math.abs(y) > 1e4) return null;
  }
  return null;
}

// ===== Klasifikasi titik kesetimbangan dari trace dan determinan =====
function cls(J){
  const tr = J[0][0] + J[1][1];
  const det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
  const D = tr*tr - 4*det, e = 1e-6;
  let t, s;
  if (Math.abs(det) < e){ t = 'Non-hiperbolik (ada λ = 0): linearisasi tidak cukup'; s = 'w'; }
  else if (det < 0){ t = 'Saddle (pelana), tidak stabil'; s = 'u'; }
  else if (D >= 0){
    if (tr < 0){ t = 'Nodal sink (node stabil)'; s = 's'; }
    else { t = 'Nodal source (node tidak stabil)'; s = 'u'; }
  }
  else if (Math.abs(tr) < e){ t = 'Center: linearisasi tidak cukup'; s = 'w'; }
  else if (tr < 0){ t = 'Spiral sink (stabil)'; s = 's'; }
  else { t = 'Spiral source (tidak stabil)'; s = 'u'; }

  const ev = D >= 0
    ? 'λ = ' + fmt((tr - Math.sqrt(D)) / 2) + ' dan ' + fmt((tr + Math.sqrt(D)) / 2)
    : 'λ = ' + fmt(tr / 2) + ' ± ' + fmt(Math.sqrt(-D) / 2) + 'i';
  return {tr, det, D, t, s, ev};
}

// Menulis sistem linearisasi, misal: u′ = π·v
const lin = (a, b, p, q) => {
  const s = [[a,p],[b,q]]
    .filter(z => z[0] !== 0)
    .map(z => (z[0] < 0 ? '− ' : '+ ') + (Math.abs(z[0]) === 1 ? '' : fmt(Math.abs(z[0])) + '·') + z[1])
    .join(' ').replace(/^\+ /, '');
  return s || '0';
};

// ===== Antarmuka =====
function pv(){
  const o = {};
  document.querySelectorAll('#params input').forEach(i => o[i.dataset.k] = +i.value);
  return o;
}

function tabs(){
  $('tabs').innerHTML = '';
  PR.forEach((p, i) => {
    const b = document.createElement('button');
    b.textContent = p.n;
    b.setAttribute('aria-pressed', i === cur);
    b.onclick = () => { cur = i; load(); };
    $('tabs').appendChild(b);
  });
}

function load(){
  const p = PR[cur];
  tabs();
  $('info').innerHTML = p.t;
  $('dim').value = p.dim;
  $('fx').value = p.f;
  $('gy').value = p.g || '';
  const d = p.d;
  $('x0').value = d[0]; $('x1').value = d[1];
  $('y0').value = d[2] ?? -1; $('y1').value = d[3] ?? 1;

  const box = $('params');
  box.innerHTML = '';
  for (const k in p.p){
    const [v, lo, hi] = p.p[k], l = document.createElement('label');
    l.innerHTML = k + ' = <span class="pv">' + v + '</span><input type="range" data-k="' + k + '" min="' + lo + '" max="' + hi + '" step="0.05" value="' + v + '">';
    l.querySelector('input').oninput = e => { l.querySelector('.pv').textContent = e.target.value; run(); };
    box.appendChild(l);
  }
  run();
}

function run(){
  $('err').textContent = '';
  try{
    dim = +$('dim').value;
    const P = pv(), pn = Object.keys(P), pvv = Object.values(P);
    $('lg').style.display = dim === 2 ? '' : 'none';
    document.querySelectorAll('.d2').forEach(e => e.style.display = dim === 2 ? '' : 'none');
    $('lf').textContent = dim === 2 ? 'f(x,y)' : 'f(x)';

    dom = [+$('x0').value, +$('x1').value, +$('y0').value, +$('y1').value];
    if (!(dom[1] > dom[0]) || (dim === 2 && !(dom[3] > dom[2])))
      throw Error('Batas maksimum harus lebih besar dari batas minimum.');

    F = comp($('fx').value, pn, pvv);
    G = dim === 2 ? comp($('gy').value, pn, pvv) : null;
    F(0,0); if (G) G(0,0);

    dim === 1 ? solve1() : solve2();
  } catch (e){
    $('err').textContent = 'Persamaan tidak valid: ' + e.message;
  }
}

// ===== Sistem 1 dimensi: x' = f(x) =====
function solve1(){
  const [a, b] = dom, f = x => F(x, 0), N = 2000, st = (b - a) / N, xs = [];
  const add = x => { if (!xs.some(z => Math.abs(z - x) < 1e-4)) xs.push(r6(x)); };

  // Akar dengan perubahan tanda (metode biseksi)
  let p = f(a);
  for (let i = 1; i <= N; i++){
    const x0 = a + (i-1)*st, x1 = a + i*st, v = f(x1);
    if (p * v < 0){
      let l = x0, h = x1;
      for (let k = 0; k < 60; k++){ const m = (l+h)/2; f(l)*f(m) <= 0 ? h = m : l = m; }
      add((l+h)/2);
    }
    p = v;
  }
  // Akar kembar (grafik menyentuh sumbu tanpa berganti tanda)
  for (let i = 1; i < N; i++){
    const x = a + i*st, q = Math.abs(f(x));
    if (q <= Math.abs(f(x-st)) && q <= Math.abs(f(x+st))){
      let l = x - st, h = x + st;
      for (let k = 0; k < 80; k++){
        const m1 = l + (h-l)/3, m2 = h - (h-l)/3;
        Math.abs(f(m1)) < Math.abs(f(m2)) ? h = m2 : l = m1;
      }
      const m = (l+h)/2;
      if (Math.abs(f(m)) < 1e-9) add(m);
    }
  }

  xs.sort((u, v) => u - v);
  eqs = xs.map(x => {
    const h = 1e-6, d = (f(x+h) - f(x-h)) / (2*h), L = f(x - .01), R = f(x + .01);
    let t, s;
    if (Math.abs(d) > 1e-4){ s = d < 0 ? 's' : 'u'; t = d < 0 ? 'Stabil' : 'Tidak stabil'; }
    else if (L > 0 && R < 0){ s = 's'; t = 'Stabil'; }
    else if (L < 0 && R > 0){ s = 'u'; t = 'Tidak stabil'; }
    else { s = 'w'; t = 'Setengah stabil (f′ = 0, linearisasi gagal)'; }
    return {x, d: r6(d), s, t};
  });

  let h = '<div class="card wide"><b>Titik tetap untuk f(x) = 0</b>';
  h += eqs.length
    ? '<table><tr><th>x*</th><th>f′(x*)</th><th>Sifat</th></tr>' +
      eqs.map(e => '<tr><td>' + fmt(e.x) + '</td><td>' + fmt(e.d) + '</td><td class="' + e.s + '"><span class="badge">' + e.t + '</span></td></tr>').join('') +
      '</table>'
    : '<p>Tidak ada titik tetap pada rentang x ini.</p>';
  $('res').innerHTML = h + '<p class="note">Aturan: f′(x*) &lt; 0 berarti stabil, f′(x*) &gt; 0 berarti tidak stabil, f′(x*) = 0 berarti linearisasi gagal.</p></div>';
}

// ===== Sistem 2 dimensi: x' = f(x,y), y' = g(x,y) =====
function solve2(){
  const [a, b, c, d] = dom, pts = [];

  // Newton-Raphson dari grid 13 x 13 tebakan awal
  for (let i = 0; i <= 12; i++) for (let j = 0; j <= 12; j++){
    const r = newton(F, G, a + (b-a)*i/12, c + (d-c)*j/12);
    if (!r) continue;
    const x = r6(r[0]), y = r6(r[1]);
    if (x < a-1e-6 || x > b+1e-6 || y < c-1e-6 || y > d+1e-6) continue;
    if (!pts.some(q => Math.hypot(q[0]-x, q[1]-y) < 1e-3)) pts.push([x, y]);
  }
  pts.sort((u, v) => u[0]-v[0] || u[1]-v[1]);

  eqs = pts.map(p => {
    const J = jac(F, G, p[0], p[1]).map(r => r.map(r6));
    return {x: p[0], y: p[1], J, ...cls(J)};
  });

  $('res').innerHTML =
    '<p class="note">Ditemukan ' + eqs.length + ' titik kesetimbangan pada rentang yang dipilih (Newton-Raphson dari 169 tebakan awal).</p>' +
    '<div class="res">' + eqs.map((e, i) => {
      const J = e.J;
      return '<div class="eq ' + e.s + '">' +
        '<h3>Titik ' + (i+1) + ': (' + fmt(e.x) + ', ' + fmt(e.y) + ')</h3>' +
        '<span class="badge">' + e.t + '</span>' +
        '<div class="jrow"><span class="jl">J =</span><div class="mat">' +
          '<span>' + fmt(J[0][0]) + '</span><span>' + fmt(J[0][1]) + '</span>' +
          '<span>' + fmt(J[1][0]) + '</span><span>' + fmt(J[1][1]) + '</span>' +
        '</div></div>' +
        '<div class="m">tr = ' + fmt(e.tr) + ', det = ' + fmt(e.det) + ', tr² − 4det = ' + fmt(e.D) + '</div>' +
        '<div class="m">' + e.ev + '</div>' +
        '<div class="m">Linearisasi: u′ = ' + lin(J[0][0], J[0][1], 'u', 'v') + ', v′ = ' + lin(J[1][0], J[1][1], 'u', 'v') + '</div>' +
      '</div>';
    }).join('') + '</div>';
}

// ===== Mulai =====
['dim','fx','gy','x0','x1','y0','y1'].forEach(i => $(i).addEventListener('input', run));
load();