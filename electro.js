// Create magnets array (global, like charges)
const magnets = [];

function electricFieldAt(x, y, charges, options = {}) {
    // options: {k:1, softening:0.1}
    const k = options.k ?? 1;
    const eps = options.softening ?? 0.1;
    let Ex = 0, Ey = 0;
    for (const c of charges) {
        const dx = x - c.x;
        const dy = y - c.y;
        const r2 = dx*dx + dy*dy + eps*eps;
        const invR3 = 1 / (Math.sqrt(r2) * r2); // 1 / r^3
        const factor = k * c.q * invR3;
        Ex += factor * dx;
        Ey += factor * dy;
    }
    return {Ex, Ey};
}

function magneticFieldAt(x, y, magnet, options = {}) {
    const soft = options.softening ?? 0.1;
    const dx = x - magnet.x;
    const dy = y - magnet.y;
    const r2 = dx*dx + dy*dy + soft*soft;
    const r = Math.sqrt(r2);
    const mx = magnet.strength * Math.cos(magnet.angle);
    const my = magnet.strength * Math.sin(magnet.angle);
    const dot = mx * dx + my * dy;
    const r5 = r2 * r2 * r;
    const r3 = r2 * r;
    const Bx = (3 * dot * dx / r5) - (mx / r3);
    const By = (3 * dot * dy / r5) - (my / r3);
    return {Bx, By};
}

function totalMagneticFieldAt(x, y, magnets, options = {}) {
    let Bx = 0, By = 0;
    for (const mag of magnets) {
        const B = magneticFieldAt(x, y, mag, options);
        Bx += B.Bx;
        By += B.By;
    }
    return {Bx, By};
}
  
function computeAccelerations(charges, magnets, options = {}) {
    const k = options.k ?? 1;
    const soft = options.softening ?? 0.1;
    const maxAccel = options.maxAccel ?? 2000;
    const N = charges.length;

    for (let i = 0; i < N; i++) {
        charges[i].ax = 0;
        charges[i].ay = 0;
    }

    for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
            const a = charges[i], b = charges[j];
            const dx = a.x - b.x;
            const dy = a.y - b.y;
            const r2 = dx*dx + dy*dy + soft*soft;
            const r = Math.sqrt(r2);
            const invR3 = 1 / (r * r2);
            const force = k * a.q * b.q * invR3;
            let fx = force * dx;
            let fy = force * dy;

            // clamp accelerations
            let ax_inc = fx / a.m;
            let ay_inc = fy / a.m;
            let magA = Math.hypot(ax_inc, ay_inc);
            if (magA > maxAccel) {
                const s = maxAccel / magA;
                ax_inc *= s;
                ay_inc *= s;
            }

            let bx_inc = -fx / b.m;
            let by_inc = -fy / b.m;
            let magB = Math.hypot(bx_inc, by_inc);
            if (magB > maxAccel) {
                const s = maxAccel / magB;
                bx_inc *= s;
                by_inc *= s;
            }

            a.ax += ax_inc;
            a.ay += ay_inc;
            b.ax += bx_inc;
            b.ay += by_inc;
        }
    }
 
    // Magnetic Lorentz force (now magnets is available)
    for (let i = 0; i < N; i++) {
        const c = charges[i];
        const B = totalMagneticFieldAt(c.x, c.y, magnets, options);
        const Bmag = Math.hypot(B.Bx, B.By);
        const vmag = Math.hypot(c.vx, c.vy);
        if (vmag > 0 && Bmag > 0) {
            const perp_x = -c.vy;
            const perp_y = c.vx;
            const perp_mag = Math.hypot(perp_x, perp_y);
            const fx = c.q * (perp_x / perp_mag) * Bmag * vmag;
            const fy = c.q * (perp_y / perp_mag) * Bmag * vmag;
            c.ax += fx / c.m;
            c.ay += fy / c.m;
        }
    }
}

function stepPhysics(charges, magnets, dt, options = {}) {
    const damping = options.damping ?? 1.0;
    const N = charges.length;

    computeAccelerations(charges, magnets, options);  // pass magnets here

    // integrate (semi-implicit Euler)
    for (let p of charges) {
        if (p.pinned) continue;
        p.vx += p.ax * dt;
        p.vy += p.ay * dt;
        p.vx *= damping;
        p.vy *= damping;
        p.x += p.vx * dt;
        p.y += p.vy * dt;
    }
}

// RK4 integrator (more accurate, conserves energy better)
function stepPhysicsRK4(charges, magnets, dt, options = {}) {
    const damping = options.damping ?? 1.0;
    const N = charges.length;

    // Save initial state
    const state0 = charges.map(c => ({
        x: c.x, y: c.y, vx: c.vx, vy: c.vy
    }));

    // k1: evaluate at t
    computeAccelerations(charges, magnets, options);  // k1
    const k1 = charges.map(c => ({
        vx: c.ax, vy: c.ay, // dv/dt = a
        x: c.vx, y: c.vy    // dx/dt = v
    }));

    // k2: evaluate at t + dt/2
    for (let i = 0; i < N; i++) {
        if (charges[i].pinned) continue;
        charges[i].x = state0[i].x + 0.5 * dt * k1[i].x;
        charges[i].y = state0[i].y + 0.5 * dt * k1[i].y;
        charges[i].vx = state0[i].vx + 0.5 * dt * k1[i].vx;
        charges[i].vy = state0[i].vy + 0.5 * dt * k1[i].vy;
    }
    computeAccelerations(charges, magnets, options);  // k2
    const k2 = charges.map(c => ({
        vx: c.ax, vy: c.ay,
        x: c.vx, y: c.vy
    }));

    // k3: evaluate at t + dt/2 again (same time, different derivative)
    for (let i = 0; i < N; i++) {
        if (charges[i].pinned) continue;
        charges[i].x = state0[i].x + 0.5 * dt * k2[i].x;
        charges[i].y = state0[i].y + 0.5 * dt * k2[i].y;
        charges[i].vx = state0[i].vx + 0.5 * dt * k2[i].vx;
        charges[i].vy = state0[i].vy + 0.5 * dt * k2[i].vy;
    }
    computeAccelerations(charges, magnets, options);  // k3
    const k3 = charges.map(c => ({
        vx: c.ax, vy: c.ay,
        x: c.vx, y: c.vy
    }));

    // k4: evaluate at t + dt
    for (let i = 0; i < N; i++) {
        if (charges[i].pinned) continue;
        charges[i].x = state0[i].x + dt * k3[i].x;
        charges[i].y = state0[i].y + dt * k3[i].y;
        charges[i].vx = state0[i].vx + dt * k3[i].vx;
        charges[i].vy = state0[i].vy + dt * k3[i].vy;
    }
    computeAccelerations(charges, magnets, options);  // k4
    const k4 = charges.map(c => ({
        vx: c.ax, vy: c.ay,
        x: c.vx, y: c.vy
    }));

    // Combine: state = state0 + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    for (let i = 0; i < N; i++) {
        if (charges[i].pinned) continue;
        const s = state0[i];
        const dxdt = (k1[i].x + 2*k2[i].x + 2*k3[i].x + k4[i].x) / 6;
        const dydt = (k1[i].y + 2*k2[i].y + 2*k3[i].y + k4[i].y) / 6;
        const dvxdt = (k1[i].vx + 2*k2[i].vx + 2*k3[i].vx + k4[i].vx) / 6;
        const dvydt = (k1[i].vy + 2*k2[i].vy + 2*k3[i].vy + k4[i].vy) / 6;

        charges[i].x = s.x + dt * dxdt;
        charges[i].y = s.y + dt * dydt;
        charges[i].vx = (s.vx + dt * dvxdt) * damping;
        charges[i].vy = (s.vy + dt * dvydt) * damping;
    }
}

// (no filepath — optional helper)
function stepPhysicsUsingFields(charges, dt, options = {}) {
    const k = options.k ?? 1;
    const soft = options.softening ?? 0.1;
    const damping = options.damping ?? 1.0;
    const N = charges.length;

    // zero
    for (let i = 0; i < N; i++) { charges[i].ax = 0; charges[i].ay = 0; }

    for (let i = 0; i < N; i++) {
        const c = charges[i];
        // compute E at c from all OTHER charges
        let Ex = 0, Ey = 0;
        for (let j = 0; j < N; j++) {
            if (i === j) continue;
            const o = charges[j];
            const dx = c.x - o.x;
            const dy = c.y - o.y;
            const r2 = dx*dx + dy*dy + soft*soft;
            const invR3 = 1 / (Math.sqrt(r2) * r2);
            const factor = k * o.q * invR3;
            Ex += factor * dx;
            Ey += factor * dy;
        }
        // F = q * E, a = F / m
        const fx = c.q * Ex;
        const fy = c.q * Ey;
        c.ax = fx / c.m;
        c.ay = fy / c.m;
    }

    // integrate (same as your code)
    for (let p of charges) {
        if (p.pinned) continue;
        p.vx += p.ax * dt;
        p.vy += p.ay * dt;
        p.vx *= damping;
        p.vy *= damping;
        p.x += p.vx * dt;
        p.y += p.vy * dt;
    }
}

// create canvas
const canvas = document.createElement('canvas');
canvas.width = 800;
canvas.height = 600;
canvas.style.display = 'block';
canvas.style.margin = '0 auto';
document.body.appendChild(canvas);
const ctx = canvas.getContext('2d');

// helper: get mouse pos in canvas coords
function getMousePos(evt) {
    const r = canvas.getBoundingClientRect();
    return { x: (evt.clientX - r.left) * (canvas.width / r.width),
             y: (evt.clientY - r.top)  * (canvas.height / r.height) };
}

// controls
const Qinput = document.getElementById('Qinput');
const Minput = document.getElementById('Minput');
const addChargeBtn = document.getElementById('addChargeBtn');
const delChargeBtn = document.getElementById('delChargeBtn');
const eulerBtn = document.getElementById('eulerBtn');
const rk4Btn = document.getElementById('rk4Btn');
const integratorLabel = document.getElementById('integratorLabel');
const addMagnetBtn = document.getElementById('addMagnetBtn');
const delMagnetBtn = document.getElementById('delMagnetBtn');
const magStrengthInput = document.getElementById('magStrengthInput');

addMagnetBtn.addEventListener('click', () => {
    const strength = parseFloat(magStrengthInput.value)*1000 || 10;
    magnets.push({
        x: canvas.width / 2,
        y: canvas.height / 2,
        angle: 0,
        strength: strength,
        pinned: false
    });    
});

delMagnetBtn.addEventListener('click', () => {
    if (magnets.length > 0) {
        magnets.pop();
    }
});

let useRK4 = false; // starts with Euler

eulerBtn.addEventListener('click', () => {
    useRK4 = false;
    integratorLabel.textContent = 'Euler';
});

rk4Btn.addEventListener('click', () => {
    useRK4 = true;
    integratorLabel.textContent = 'RK4';
});

addChargeBtn.addEventListener('click', () => {
    const q = parseFloat(Qinput.value) || 1;
    const m = parseFloat(Minput.value) || 1;

    const newCharge = { 
        x: canvas.width / 2,   // center x
        y: canvas.height / 2,  // center y
        vx: 0, vy: 0, q: q, m: m, pinned: false 
    };
    charges.push(newCharge);
});

delChargeBtn.addEventListener('click', () => {
    if (charges.length > 0) {
        charges.pop();
    }
});

// draw function: optional field arrows + charges
function draw(charges, magnets, options = {}) {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    // background
    ctx.fillStyle = '#101010';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // draw coarse E-field arrows (optional)
    if (options.drawField) {
        const spacing = options.gridSpacing || 32;
        ctx.strokeStyle = 'rgba(255,255,255,0.6)';
        ctx.lineWidth = 1;
        for (let y = spacing/2; y < canvas.height; y += spacing) {
            for (let x = spacing/2; x < canvas.width; x += spacing) {
                const E = electricFieldAt(x, y, charges, {k: options.k, softening: options.softening});
                const mag = Math.hypot(E.Ex, E.Ey);
                if (mag === 0) continue;
                // scale arrow for visibility
                const s = options.fieldScale ?? 20;
                const ex = (E.Ex / mag) * Math.min(mag * s, spacing*0.9);
                const ey = (E.Ey / mag) * Math.min(mag * s, spacing*0.9);
                drawArrow(ctx, x, y, x + ex, y + ey);
            }
        }
    }

    // draw charges
    for (const c of charges) {
        const r = c._r ?? Math.max(6, Math.min(20, Math.abs(c.q) * 4));
        ctx.beginPath();
        ctx.fillStyle = c.q > 0 ? '#ff4444' : '#4444ff';
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1.5;
        ctx.arc(c.x, c.y, r, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
        // label
        ctx.fillStyle = '#fff';
        ctx.font = '12px sans-serif';
        ctx.fillText((c.q>0?'+':'')+c.q.toString(), c.x + r + 4, c.y + 4);
    }

    // magnets with colored poles
    for (const m of magnets) {
        const w = 40, h = 20;
        ctx.save();
        ctx.translate(m.x, m.y);
        ctx.rotate(m.angle);
        
        // North pole (left, red)
        ctx.fillStyle = '#ff6666';
        ctx.fillRect(-w/2, -h/2, w/2, h);
        
        // South pole (right, blue)
        ctx.fillStyle = '#6666ff';
        ctx.fillRect(0, -h/2, w/2, h);
        
        // outline
        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 2;
        ctx.strokeRect(-w/2, -h/2, w, h);
        
        // labels
        ctx.fillStyle = '#fff';
        ctx.font = 'bold 12px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('N', -w/4, 4);
        ctx.fillText('S', w/4, 4);
        
        ctx.restore();
    }
}

// small helper that draws an arrow from (x1,y1) to (x2,y2)
function drawArrow(ctx, x1, y1, x2, y2) {
    const dx = x2 - x1, dy = y2 - y1;
    const ang = Math.atan2(dy, dx);
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2, y2);
    ctx.stroke();
    // arrow head
    const ah = 6;
    ctx.beginPath();
    ctx.moveTo(x2, y2);
    ctx.lineTo(x2 - ah * Math.cos(ang - 0.4), y2 - ah * Math.sin(ang - 0.4));
    ctx.lineTo(x2 - ah * Math.cos(ang + 0.4), y2 - ah * Math.sin(ang + 0.4));
    ctx.closePath();
    ctx.fillStyle = ctx.strokeStyle;
    ctx.fill();
}

// simple dragging logic
let dragging = null;
let draggingType = null;

canvas.addEventListener('mousedown', (e) => {
    const pos = getMousePos(e);

    // pick nearest charge within radius
    let best = null, bestDist = Infinity;
    for (const c of charges) {
        const r = c._r ?? Math.max(6, Math.min(20, Math.abs(c.q) * 4));
        const d2 = (c.x - pos.x)**2 + (c.y - pos.y)**2;
        if (d2 < (r + 6)**2 && d2 < bestDist) { best = c; bestDist = d2; }
    }

    for (const m of magnets) {
        const w = 40, h = 20;
        const d2 = (m.x - pos.x)**2 + (m.y - pos.y)**2;
        if (d2 < (Math.max(w,h)/2 + 6)**2 && d2 < bestDist) { best = m; bestDist = d2; draggingType = 'magnet'; } 
    }

    if (best) {
        dragging = best;
        dragging.pinned = true;
        dragging.vx = dragging.vy = 0;
    } else {
        // optional: add new charge on click (example)
        // charges.push({x:pos.x, y:pos.y, vx:0, vy:0, q:1, m:1, pinned:false});
    }
});

canvas.addEventListener('mousemove', (e) => {
    if (!dragging) return;
    const pos = getMousePos(e);
    dragging.x = pos.x;
    dragging.y = pos.y;
});

canvas.addEventListener('mouseup', () => { 
    if (dragging) dragging.pinned = false; 
    dragging = null;
    draggingType = null;
});

canvas.addEventListener('mouseleave', () => { 
    if (dragging) dragging.pinned = false; 
    dragging = null;
    draggingType = null;
});

// animation loop using your stepPhysics
let last = performance.now();
// increase k so accelerations are visible; lower if things blow up
const simOptions = {k: 10000, softening:0.4, damping:0.999, maxAccel:2000};
const drawOptions = {drawField:true, gridSpacing:28, fieldScale:18, k:1, softening:0.4};

function loop(now) {
    const dt = Math.min(0.03, (now - last) / 1000);
    last = now;
    
    if (useRK4) {
        stepPhysicsRK4(charges, magnets, dt, simOptions);  // add magnets
    } else {
        stepPhysics(charges, magnets, dt, simOptions);  // add magnets
    }
    
    draw(charges, magnets, drawOptions);  // add magnets
    requestAnimationFrame(loop);
}
requestAnimationFrame(loop);