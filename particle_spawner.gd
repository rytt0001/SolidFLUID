extends Node2D

@export 
var maxNumOfParticlesToSpawn := 0;
var curNumOfParticles:= 0;

@export 
var SpawnPerSecond := 10.0;

var curTimer := 0.0;

@export
var ParticleInstance :=preload("res://particle_instance_body.tscn");
# Called when the node enters the scene tree for the first time.
func _ready():
	pass # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta):
	curTimer+=delta;
	if(curTimer >= 1.0/SpawnPerSecond && curNumOfParticles < maxNumOfParticlesToSpawn):
		#spawn
		var new_object = ParticleInstance.instantiate();
		get_tree().root.add_child(new_object)
		new_object.position = position
		new_object.position.x += randf_range(-10.0,10)
		curNumOfParticles= curNumOfParticles+1;
		curTimer = 0.0;
		
	pass
