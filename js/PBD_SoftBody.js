
/*
Shape Matching for Position Based Dynamics
based on the paper from Matthias Mueller:
 http://matthias-mueller-fischer.ch/talks/2017-EG-CourseNotes.pdf
 */


import {Particle} from './Particle.js';

const radius = 6; //Visual radius of the particles
const deltaTime = 0.06; //speed of the simulation
const iterations = 30; //Number of iterations for the shape matching solver
const gravity = 10;


let shapes = [];
let particles = [];
let ii = 0;
let container = document.querySelector('.container');
let looping = true;


function evaluateCenterOfMass(particles) {
    let centerOfMass = {x: 0, y:0}
    for(let particle of particles) {
        centerOfMass.x += particle.position.x;
        centerOfMass.y += particle.position.y;
    }
    centerOfMass.x /= particles.length;
    centerOfMass.y /= particles.length;
    return centerOfMass;
}


/*
The name in the shape and the particles is used to allocate the particles in the DOM (to place them in their positions),
but most important it's used for collisions, where a single pool of particles will be used to define the neighborhood
search for close collisions, the particles from the same shape shouldn't be evaluated for collisions (avoiding inner
collisions).

cx... center X
cy... center Y
stiffness... controls the soft body softness, the range goes from [0 - 1], values above 0.4 are advised.
 */

let shapeDefinition = [
    {radius: 56, cx: window.innerWidth * 0.5, cy: window.innerHeight * 0.5, stiffness: 0.005}
]


//create the particles and shapes.
for(let shapeDef of shapeDefinition) {
    particles = [];

    //actual shape, can be a regular mesh.
    let velocityDirection = 2 * Math.PI * Math.random();
    for(let i = 0; i < 60; i ++) {
        const angle = 2 * Math.PI * i / 60;
        const r = shapeDef.radius;
        particles.push(new Particle(`shape${ii}`, document.body, r * Math.cos(angle) + shapeDef.cx, r * Math.sin(angle) + shapeDef.cy, particles.length, velocityDirection));
        particles[particles.length - 1].radius(radius);
    }

    shapes[ii] = {radius: Math.sqrt(90), name: `shape${ii}`, particles: [...particles], initialCenterOfMass: evaluateCenterOfMass(particles), stiffness: shapeDef.stiffness, centerOfMass: evaluateCenterOfMass(particles)}

    ii++;
}

//TODO this is precalculated should be saved in a texture. Should check if it's better to generate this in the CPU or the GPU
//Evaluate the relative position for each particle relative to the center of mass to each shape
for(let shape of shapes) {

    //Evaluate Aqq for the linear extension
    let Aqq = mat2.fromValues(0, 0, 0, 0);

    for(let particle of shape.particles) {
        particle.initialRelativePosition.x = particle.position.x - shape.initialCenterOfMass.x;
        particle.initialRelativePosition.y = particle.position.y - shape.initialCenterOfMass.y;

        const q = particle.initialRelativePosition;
        let pqi = mat2.fromValues(q.x * q.x, q.y * q.x, q.x * q.y, q.y * q.y);

        mat2.add(Aqq, Aqq, pqi);
    }

    mat2.invert(Aqq, Aqq);
    shape.Aqq = Aqq;
}

particles = null;



/*
 Collision evaluation between particles and borders
 This is done particle wise and it could be O(n2) without
 any acceleration:
 - Acceleration is done evaluating bounding volumes / shapes
 - For the GPU a neighborhood search could be implemented  <<--- implies Harada sorting in buckets.
 */

let evaluateCollisions = (shape) => {

    container.style.borderBottom = "solid 3px #0000ff";
    container.style.borderTop = "solid 3px #0000ff";
    container.style.borderLeft = "solid 3px #0000ff";
    container.style.borderRight = "solid 3px #0000ff";

    for(const particle of shape.particles) {


        const center = {x:window.innerWidth * 0.495, y: window.innerHeight * 0.495}; //center of the border box

        //Collision on box
        const boxSize = {x: window.innerWidth * 0.495, y: window.innerHeight * 0.495}

        const xLocal = {x: particle.position.x - center.x, y: particle.position.y - center.y}

        const contactLocalPoint = {x: Math.min(Math.max(xLocal.x, - boxSize.x), boxSize.x),
                                 y: Math.min(Math.max(xLocal.y, - boxSize.y), boxSize.y)}

        const contactPoint = {x: contactLocalPoint.x + center.x, y: contactLocalPoint.y + center.y}

        const vector = {x: contactPoint.x - particle.position.x, y: contactPoint.y - particle.position.y}
        const dist = Math.sqrt(vector.x * vector.x + vector.y * vector.y);

        if(dist > 0.) {

            particle.position.x = contactPoint.x;
            particle.position.y = contactPoint.y;

            let normal = {x: 0, y: 0}

            //Touches the left wall
            if(particle.position.x == 0 && particle.position.y > 0) {
                normal.x = 1;
                container.style.borderLeft = "solid 3px #ff0000";
            }

            //Touches the right wall
            if(particle.position.x == 2 * window.innerWidth * 0.495 && particle.position.y > 0) {
                normal.x = -1;
                container.style.borderRight = "solid 3px #ff0000";
            }

            //Touches the lower wall
            if(particle.position.x > 0 && particle.position.y == 0) {
                normal.y = -1;
                container.style.borderTop = "solid 3px #ff0000";
            }

            //Touched the upper wall
            if(particle.position.x > 0 && particle.position.y == 2 * window.innerHeight * 0.495) {
                normal.y = 1;
                container.style.borderBottom = "solid 3px #ff0000";
            }

        }

    }
}



/*
This function is used to evaluate the rotation matrix used to update the particles
it's called per shape, but requires to sum all the particles.
 */
let generateShapeMatrix = shape => {

    //Evaluate the relative position for each particle relative to the center of mass O(n)
    let Apq = mat2.fromValues(0, 0, 0, 0);
    for(let particle of shape.particles) {
        particle.relativePosition.x = particle.position.x - shape.centerOfMass.x;
        particle.relativePosition.y = particle.position.y - shape.centerOfMass.y;

        const p = particle.relativePosition;
        const q = particle.initialRelativePosition;
        let pqi = mat2.fromValues(p.x * q.x, p.y * q.x, p.x * q.y, p.y * q.y);
        mat2.add(Apq, Apq, pqi);
    }

    //Calculate the transpose of Apq
    let ApqT = mat2.create();
    mat2.transpose(ApqT, Apq);

    //Calculate S
    let S = mat2.create();
    mat2.mul(S, ApqT, Apq);

    /*
     Using Babylonian method to evaluate the square root.
     */

    let Si = mat2.create();
    let SiI = mat2.create();
    for(let i = 0; i < 50; i ++) {
        SiI = mat2.invert(SiI, Si);
        mat2.multiplyScalar(Si, mat2.add(Si, Si, mat2.multiply(SiI, S, SiI)), 0.5);
    }

    let Sinv = mat2.create();
    mat2.invert(Sinv, Si);

    //Calculate rotation matrix
    let R = mat2.create();
    mat2.mul(R, Apq, Sinv);


    //Linear extension
    let A = mat2.create();
    A = mat2.mul(A, Apq, shape.Aqq);
    let detA = mat2.determinant(A);
    detA = Math.pow(detA, 1/ 3);
    mat2.multiplyScalar(A, A, detA);

    const beta = 0.2;
    mat2.add(A, mat2.multiplyScalar(A, A, beta), mat2.multiplyScalar(R, R, 1 - beta));

    return A;
}



//Calculate the displacement for all the particles of the shape
let updateParticlesIterationPositions = (shape, A) => {
    for (let particle of shape.particles) {
        const d = particle.initialRelativePosition;
        let g = vec2.fromValues(d.x, d.y);

        //Apply the rotation for the initial relative position.
        vec2.transformMat2(g, g, A);

        particle.position.x += shape.stiffness * (g[0] + shape.centerOfMass.x - particle.position.x);
        particle.position.y += shape.stiffness * (g[1] + shape.centerOfMass.y - particle.position.y);

    }
}



let iterativeStep = (shape) => {

    for(let shape of shapes) shape.centerOfMass = evaluateCenterOfMass(shape.particles); //<<<<<<<--------- O(n)


    evaluateCollisions(shape); //<<<<---------- O(n2)


    const A = generateShapeMatrix(shape); //<<<<<<<<<--------- O(n)


    updateParticlesIterationPositions(shape, A);
}

let simulationStep = () => {

    if(looping) requestAnimationFrame(simulationStep);

    for (let shape of shapes) {

        for(let particle of shape.particles) {
//            particle.velocity.y += deltaTime * gravity;
            particle.position.x += particle.velocity.x * deltaTime;
            particle.position.y += particle.velocity.y * deltaTime;
        }

        //Solve the iterative position constrain
        for(let i = 0; i < iterations; i ++) iterativeStep(shape);


        for(let particle of shape.particles) {

            particle.x = particle.position.x;
            particle.y = particle.position.y;

            particle.velocity.x = (particle.x - particle.prevPosition.x) / deltaTime;
            particle.velocity.y = (particle.y - particle.prevPosition.y) / deltaTime;

            particle.update(particle.x, particle.y);
        }
    }
}

simulationStep();

document.body.addEventListener("click", () => {
        looping = !looping;
        if(looping) simulationStep();
    }
);

document.body.addEventListener('keypress', (e) => {
    if(e.key == "n" && !looping) simulationStep();
});
