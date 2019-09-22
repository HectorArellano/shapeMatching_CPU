

class Particle {

    constructor(name, container, x, y, id) {
        this.constrain = 0;
        this.density = 0;
        this.lambda = 0;
        this.mass = 0.1;

        this.initialRelativePosition = {x: 0, y: 0}
        this.relativePosition = {x: 0, y: 0}

        this.x = x;
        this.y = y;

        this.position = {x:this.x, y:this.y};
        this.prevPosition = {x:this.x, y:this.y};

        this.displacement = {x:0, y:0};
        this.velocity = {x:0, y:0}
        this.colorSet = false;
        this.name = name;
        this.container = container;
        this.baseId = id;

        this.id = "particle_" + String(name + "_" + id);
        let graphic = document.createElement('div');
        graphic.id = "particle_" + String(name + "_" + id);
        graphic.style.border = "#ff0000 1px solid";
        graphic.style.position = "absolute";
//        graphic.style.cursor = "pointer";
//        graphic.style.backgroundColor = "#ff0000";


//        graphic.addEventListener("mouseover", () => {
//            console.log("=================================");
//            console.log("id: " + this.id);
//            console.log("constrain: " + this.constrain);
//            console.log("lambda: " + this.lambda);
//            console.log("displacement: x=" + this.displacement.x + " y= " + this.displacement.y);
//            console.log("velocity: x=" + this.velocity.x + " y= " + this.velocity.y);
//            console.log("---------------------------------")
//
//        });

        container.appendChild(graphic);

        this.radius(4);
        this.update(this.x, this.y);
    }

    update(x, y) {

        this.position.x = x;
        this.position.y = y;

        this.prevPosition.x = x;
        this.prevPosition.y = y;

        this.x = x;
        this.y = y;

        let graphic = document.querySelector("#" + this.id);
        graphic.style.top = String(y * 8) + "px";
        graphic.style.left = String(x * 8) + "px"
    }

    distance(particle) {
        let x = this.x - particle.x;
        let y = this.y - particle.y;
        return Math.sqrt(x * x + y * y);
    }

    radius(m) {
        let r = m - 1;
        let graphic = document.querySelector("#" + this.id);
        graphic.style.width = String(r) + "px";
        graphic.style.height = String(r) + "px";
        graphic.style.borderRadius = String(r) + "px";
    }

}

export {Particle}