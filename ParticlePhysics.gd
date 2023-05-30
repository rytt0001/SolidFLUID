extends RigidBody2D

var linear_acceleration  := Vector2(0,0);
var gravity := Vector2.DOWN * 9.80
var coef_mu := 0.2 # frottement
var coef_e := 0.9 # restitution
@export
var do_collision_maths := true
var selfScript 


var lnormal := Vector2(0,0)
var lobjnormal :=Vector2(0,0)
# Called when the node enters the scene tree for the first time.
func _ready():
	add_to_group("Particle")
	pass # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta):
	
	pass

func _draw():
	draw_line(Vector2(0,0), Vector2(0,0)+(lnormal*15), Color.RED, 2)
	#draw_line(Vector2(0,0), Vector2(0,0)+(lobjnormal*15), Color.RED, 2)
	
	pass
func _physics_process(delta):
	
	#return
	
	custom_apply_force(gravity*gravity_scale);
	linear_velocity += linear_acceleration * delta;
	
	
	
	var result =move_and_collide(linear_velocity*delta, true);
	get_colliding_bodies()
	if (result != null && do_collision_maths):
		var obj : RigidBody2D
		obj =result.get_collider()
		var objNormal = -result.get_normal()
		lnormal=result.get_depth()*objNormal
		queue_redraw()
		var m1 = mass
		var v1 = linear_velocity
		var v2 = obj.linear_velocity
		var m2 = obj.mass
		var vRel = (v1-v2).dot(result.get_normal())
		
		
		#POS Correction
		var damping = 0.7;
		var correction = (result.get_depth() * damping)/((1/m1) + (1/m2));
		position += (1/m1) * correction * result.get_normal();
		if (obj.is_in_group("Particle")):
			obj.position -= (1/m2) * correction * result.get_normal();
		
		if(vRel < 0.0):
			
			var J = (-(1.0+coef_e)*vRel)/(1.0/m1 + 1.0/m2)
			
			var colImpulse1 =J * (1.0/m1) * result.get_normal();
			var colImpulse2 = J * (1.0/m2) * result.get_normal();
			var tangent = result.get_normal().rotated(deg_to_rad(90))
			tangent = tangent.normalized()
			var vtan = (v1-v2).dot(tangent)
			var J2 = -vtan/((1.0/m1)+(1.0/m2));
			var JA2 = clampf(J2,-absf(J)*coef_mu, absf(J)*coef_mu);
			
			
			
			linear_velocity += colImpulse1
			
			linear_velocity+= JA2 *tangent * (1.0/m1)
			if(obj.is_in_group("Particle")):
				#obj.do_collision_maths = false
				var JB2 = clampf(J2,-absf(J)*obj.coef_mu, absf(J)*obj.coef_mu);
				obj.linear_velocity -= colImpulse2
				obj.linear_velocity -= JB2 *tangent * (1.0/m2)
		
		
		
		
			if (obj.has_meta("is_wall")):
				obj.linear_velocity = Vector2(0,0)
	else:
		do_collision_maths = true;
	
	position += linear_velocity*delta;
	
	linear_acceleration = Vector2(0,0);
	pass
func custom_apply_force(force : Vector2):
	linear_acceleration += force / mass;
	pass
	

