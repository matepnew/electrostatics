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

function stepPhysics(charges, dt, options = {}) {
    const k = options.k ?? 1;
    const soft = options.softening ?? 0.1;
    const damping = options.damping ?? 1.0; // e.g., 0.999
    const maxAccel = options.maxAccel ?? 2000; // safety clamp to avoid explosions

    // compute forces (naive O(N^2))
    const N = charges.length;
    // zero accelerations
    for (let i = 0; i < N; i++) {
        charges[i].ax = 0;
        charges[i].ay = 0;
    }
    for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
            const a = charges[i], b = charges[j];
            // Use vector r_a - r_b so sign of force (attract/repel) is correct:
            const dx = a.x - b.x;
            const dy = a.y - b.y;
            const r2 = dx*dx + dy*dy + soft*soft;
            const r = Math.sqrt(r2);
            const invR3 = 1 / (r * r2); // 1 / r^3
            const force = k * a.q * b.q * invR3;
            // force vector on a due to b
            let fx = force * dx;
            let fy = force * dy;

            // convert to accelerations and clamp magnitude for stability
            let ax_inc = fx / a.m;
            let ay_inc = fy / a.m;
            let magA = Math.hypot(ax_inc, ay_inc);
            if (magA > maxAccel) {
                const s = maxAccel / magA;
                ax_inc *= s;
                ay_inc *= s;
            }

            let bx_inc = -fx / b.m; // equal and opposite on b
            let by_inc = -fy / b.m;
            let magB = Math.hypot(bx_inc, by_inc);
            if (magB > maxAccel) {
                const s = maxAccel / magB;
                bx_inc *= s;
                by_inc *= s;
            }

            // update accelerations: a += F/m
            a.ax += ax_inc;
            a.ay += ay_inc;
            b.ax += bx_inc;
            b.ay += by_inc;
        }
    }
    // integrate
    for (let p of charges) {
        if (p.pinned) continue; // being dragged by user
        p.vx += p.ax * dt;
        p.vy += p.ay * dt;
        p.vx *= damping;
        p.vy *= damping;
        p.x += p.vx * dt;
        p.y += p.vy * dt;
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
function draw(charges, options = {}) {
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
canvas.addEventListener('mousedown', (e) => {
    const pos = getMousePos(e);
    // pick nearest charge within radius
    let best = null, bestDist = Infinity;
    for (const c of charges) {
        const r = c._r ?? Math.max(6, Math.min(20, Math.abs(c.q) * 4));
        const d2 = (c.x - pos.x)**2 + (c.y - pos.y)**2;
        if (d2 < (r + 6)**2 && d2 < bestDist) { best = c; bestDist = d2; }
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
canvas.addEventListener('mouseup', () => { if (dragging) dragging.pinned = false; dragging = null; });
canvas.addEventListener('mouseleave', () => { if (dragging) dragging.pinned = false; dragging = null; });

// animation loop using your stepPhysics
let last = performance.now();
// increase k so accelerations are visible; lower if things blow up
const simOptions = {k: 10000, softening:0.4, damping:0.999, maxAccel:2000};
const drawOptions = {drawField:true, gridSpacing:28, fieldScale:18, k:1, softening:0.4};

function loop(now) {
    const dt = Math.min(0.03, (now - last) / 1000); // clamp dt for stability
    last = now;
    stepPhysics(charges, dt, simOptions); // uses your function
    draw(charges, drawOptions);
    requestAnimationFrame(loop);
}
requestAnimationFrame(loop);