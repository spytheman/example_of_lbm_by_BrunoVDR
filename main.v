/**
	This program is about D2Q9 Lattice Boltzmann Method:
	See : https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods

	It's a pet project in order to use V language: https://vlang.io/

	The simulation is single threaded, probably buggy and should not be used
    for serious things.  It's very sensible to tau parameter that should be
    carefully set. This parameter is related to fluid viscosity, and is set so
    that the fluid speed doesn't exceed a speed limit that breaks simulation.
    Too narrow passage (speed increased) may reach this limit.

    profiles files MUST be of the same size of defined width and weight and
    should be 8bits per pixels. Every non zero value is considered as an
    obstacle.

    to compile the program from within source directory:

    v -prod .

	or if you want gcc as compiler:

    v -prod -cc gcc .

    This program is released under MIT license.
 */
module main

import gg
import os
import stbi
import time

const tau = 0.6 // Relaxation time, related to fluid viscosity.
const rho0 = 32.0 // Average normalised density.
const width = 512
const height = 128
const len_in_bytes = width * height * sizeof(u32)
const obstacle_color = u32(0xFF004000)
const low_color = u32(0xFFFFFF00)
const middle_color = u32(0xFF000000)
const high_color = u32(0xFF00FFFF)

// Type for rendering methods, used as parameter.
type Renderer = fn (l Lattice, cm []u32, mut output []u32)

const renderers = [vorticity, v_speed, h_speed, densities]
const renderers_names = ['vorticity', 'vertical speed', 'horizontal speed', 'density']

@[heap]
struct App {
mut:
	gg            &gg.Context = unsafe { nil }
	iidx          int
	render_method Renderer = unsafe { nil }
	pixel_buffer  []u32
	cm            []u32
	src           Lattice
	dst           Lattice
}

fn (mut app App) on_init() {
	app.iidx = app.gg.new_streaming_image(width, height, 4, pixel_format: .rgba8)
}

fn (mut app App) on_keydown(code gg.KeyCode, mod gg.Modifier, data voidptr) {
	match code {
		.v {
			// Show next view
			mut i := renderers.index(app.render_method)
			i = (i + 1) % renderers.len
			app.render_method = renderers[i]
			println('Rendering : ${renderers_names[i]}')
		}
		.escape {
			app.gg.quit()
		}
		else {}
	}
}

fn (mut app App) on_frame() {
	mut stop_watch := time.new_stopwatch()
	app.src.move(mut app.dst)
	app.dst.collide()
	app.render_method(app.dst, app.cm, mut app.pixel_buffer) // render_method can point different method !
	draw_colormap(app.cm, mut app.pixel_buffer)

	// swap src and dst buffers.
	tmp := app.src
	app.src = app.dst
	app.dst = tmp

	mut istream_image := app.gg.get_cached_image_by_idx(app.iidx)
	istream_image.update_pixel_data(unsafe { &u8(app.pixel_buffer.data) })
	stop_watch.stop()
	// println('Frame ${app.gg.frame}, loop : ${stop_watch.elapsed().milliseconds()} milliseconds. ')

	size := gg.window_size()
	app.gg.begin()
	app.gg.draw_image(0, 0, size.width, size.height, istream_image)
	app.gg.end()
}

fn main() {
	argv := os.args.len
	if argv != 2 {
		println('Usage: lbm profile_file.png')
		println('       e.g:  ./lbm profiles/circle.png')
		println('During simulation press "v" to show different parameters.')
		return
	}

	profile := stbi.load(os.args[1])!
	// Check size compatibility.
	if profile.width != width || profile.height != height {
		eprintln('Error, "${os.args[1]}" profile image must match lbm lattice size : ${profile.width}x${profile.height}')
		return
	}
	mut profile_pixels := []u8{len: profile.width * profile.height}
	for i in 0 .. profile_pixels.len {
		profile_pixels[i] = unsafe { (&u8(profile.data))[4 * i] }
	}

	mut app := &App{}
	app.gg = gg.new_context(
		width:        width * 2
		height:       height * 2
		window_title: 'Lattice Boltzmann Method [D2Q9]'
		init_fn:      app.on_init
		frame_fn:     app.on_frame
		keydown_fn:   app.on_keydown
		user_data:    app
	)

	// Build a colormap to be used
	app.cm = Colormap.dual(low_color, middle_color, high_color, 384)

	// Now create Lattices, with respect to loaded profile.
	app.src = Lattice.new(width, height, profile_pixels.data)

	app.src.add_flow(1.0, Vi.east)
	app.src.randomize(0.2)
	app.src.normalize()

	app.dst = Lattice.new(width, height, profile_pixels.data)

	// Allocate pixel buffer to draw in.
	app.pixel_buffer = []u32{len: width * height} // Dyn array heap allocated.
	app.render_method = vorticity
	println('Showing vorticiy. Press "v" to show different parameters.')
	app.gg.run()
}

fn draw_colormap(cm []u32, mut data []u32) {
	data[0] = 0xFF000000
	data[width] = 0xFF000000
	data[2 * width] = 0xFF000000
	for i in 0 .. cm.len {
		data[i + 1] = 0xFF000000
		data[i + width + 1] = cm[i]
		data[i + 1 + 2 * width] = 0xFF000000
	}
	data[cm.len] = 0xFF000000
	data[width + cm.len] = 0xFF000000
	data[2 * width + cm.len] = 0xFF000000
}

// densities is a Renderer type function
fn densities(l Lattice, cm []u32, mut output []u32) {
	mut ind := 0

	min_rho := l.min_rho()
	max_rho := l.max_rho()
	linear := (max_rho - min_rho) / (cm.len - 1)

	for c in l.m {
		if c.obstacle == true {
			output[ind] = obstacle_color
			ind++
			continue
		}

		rho := int((c.rho() - min_rho) / linear)
		output[ind] = cm[rho]
		ind++
	}
}

// h_speed is a Renderer type function
fn h_speed(l Lattice, cm []u32, mut output []u32) {
	mut ind := 0

	min_ux := l.min_ux()
	max_ux := l.max_ux()
	linear := (max_ux - min_ux) / (cm.len - 1)

	for c in l.m {
		if c.obstacle == true {
			output[ind] = obstacle_color
			ind++
			continue
		}

		rho := int((c.ux() - min_ux) / linear)
		output[ind] = cm[rho]
		ind++
	}
}

// h_speed is a Renderer type function
fn v_speed(l Lattice, cm []u32, mut output []u32) {
	mut ind := 0

	min_uy := l.min_uy()
	max_uy := l.max_uy()
	linear := (max_uy - min_uy) / (cm.len - 1)

	for c in l.m {
		if c.obstacle == true {
			output[ind] = obstacle_color
			ind++
			continue
		}

		rho := int((c.uy() - min_uy) / linear)
		output[ind] = cm[rho]
		ind++
	}
}

// vorticity is a Renderer type function
fn vorticity(l Lattice, cm []u32, mut output []u32) {
	mut min := 0.0
	mut max := 0.0
	cm2 := u32(cm.len / 2)
	mut vorticity_table := []f64{len: l.w * l.h}

	for y in 0 .. l.h {
		for x in 0 .. l.w {
			out := (y * l.w) + x
			if l.m[out].obstacle {
				vorticity_table[out] = -1000000.0
			} else {
				v := l.vorticity(x, y)
				vorticity_table[out] = v
				if min > v {
					min = v
				}
				if max < v {
					max = v
				}
			}
		}
	}

	linear := (max - min) / f64(cm.len - 1)

	for ind, v in vorticity_table {
		if v < -100.0 {
			output[ind] = obstacle_color
		} else {
			mut id := cm2 + u32(v / linear)

			if id < 0 {
				id = 0
			} else if id >= cm.len {
				id = u32(cm.len - 1)
			}
			output[ind] = cm[id]
		}
	}
}
