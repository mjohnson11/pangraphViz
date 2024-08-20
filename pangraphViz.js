
/* TODO 
- make gaps in core genome for junctions using local pangenome size
- make an alignment viewer for when you click on a block
- work on annotation
- add original e coli dataset
- general clean up
*/

const path_view_frac = 0;//0.25
const vert_bez = 0.05; //vertica bezier offset
const core_only = true;
const remove_dups = false;
const draw_core_only = false;
const show_paths_always = true;
const hover_stroke_width = core_only ? 4 : 2;
const colors = {
  'forward_path': 'black',
  'rev_path': 'black'
}

const seaborn_colorblind = [
  "#0173b2",  //A strong blue
  "#de8f05",  //A muted orange
  "#029e73",  //A greenish teal
  "#d55e00",  //A vivid red-orange
  "#cc78bc",  //A soft purple
  "#ca9161",  //A light brown or tan
  "#fbafe4",  //A pinkish color
  "#949494",  //A neutral grey
  "#ece133",  //A bright yellow
  "#56b4e9"   //A light blue
]

const greys = [
  '#AAA',
  '#666'
]

let focal_strain;
let second_strain;
let circlegraph;
let graph_data;
let treeObj;
let clickedJunction;

const dset_names = [
  'kpneumoniae_ST512',
  'mtuberculosis_ST276',
  'paeruginosa_ST235',
  'ecoli_ST69',
  'saureus_ST1'
]

const dsets = dset_names.reduce((td, n) => {
  td[n] = {
    'treefile': './example_data/'+n+'/coretree.nwk',
    'graphfile': './example_data/annotated_'+n+'.json'
  }
  return td;
}, {})

dsets['pneumo'] = {
  'treefile': '../spneumo/spneumoniae_gpsc1_clean.nwk',
  'graphfile': '../spneumo/polished_spneumo_GPSC1_graph_annot.json'
}

const use = 'saureus_ST1';
const use_test_data = false;
let annotated = true;


function show_about() {
  var about_t = d3.select("#about_text");

  // Toggle the display property between 'none' and 'block'
  if (about_t.style("display") === "block") {
      about_t.style("display", "none");
  } else {
      about_t.style("display", "block");
  }

  // Scroll the page to the #about_text element
  var about_t_position = about_t.node().getBoundingClientRect().top + window.scrollY;
  
  d3.select("html, body")
    .transition()
    .duration(1000)
    .tween("scroll", function() {
        var i = d3.interpolateNumber(window.scrollY, about_t_position);
        return function(t) { window.scrollTo(0, i(t)); };
    });
}

function point_on_circle(cx, cy, r, theta) {
  return [cx+Math.cos(theta)*r, cy+Math.sin(theta)*r];
}

function extend_off_circle(cx, cy, p, extend) {
  return [(p[0]-cx)*extend+cx, (p[1]-cy)*extend+cy];
}

let tooltip;

function setup_tooltip() {
  tooltip = d3.select('#svg_div').append('div')
    .attr('class', 'WOW_tooltip')
    .html('<h2>yeah</h2><p>uhhuh</p>');
}

function show_tooltip(x, y, text) {
  tooltip.style('left', String(x+20)+'px');
  tooltip.style('display', 'block').html(text);
  tooltip.style('top', String(y+10)+'px');
}

function hide_tooltip() {
  tooltip.style('display', 'none');
}

// Newick Parser Function
function parseNewick(newick) {
  var ancestors = [];
  var tree = {};
  var tokens = newick.split(/\s*(;|\(|\)|,|:)\s*/);
  for (var i = 0; i < tokens.length; i++) {
      var token = tokens[i];
      switch (token) {
          case '(': // new branchset
              var subtree = {};
              tree.children = [subtree];
              ancestors.push(tree);
              tree = subtree;
              break;
          case ',': // another branch
              var subtree = {};
              ancestors[ancestors.length - 1].children.push(subtree);
              tree = subtree;
              break;
          case ')': // end current branchset
              tree = ancestors.pop();
              break;
          case ':': // ignore lengths
              break;
          default:
              var x = tokens[i - 1];
              if (x == ')' || x == '(' || x == ',') {
                  tree.name = token;
              } else if (x == ':') {
                  tree.length = parseFloat(token);
              }
      }
  }
  return tree;
}

function max_tree_length(tree) {
  let tlen = tree.length || 0; // case for outer tree
  if (tree.children) {
    return tlen + Math.max(...tree.children.map( max_tree_length));
  } else {
    return tree.length;
  }
}

function tree_nodes(tree) {
  if (tree.children) {
    tree.node_count = tree.children.map(tree_nodes).reduce((summer, count) => {
      summer += count;
      return summer;
    }, 0);
  } else {
    tree.node_count = 1
  }
  return tree.node_count;
}

class SimpleTree {

  constructor(newick_string, w, h, right_buffer) {
    this.newick_string = newick_string;
    this.tree = parseNewick(this.newick_string);
    this.tree_len = max_tree_length(this.tree);
    this.tree_nodes = tree_nodes(this.tree);
    this.strain_locations = {};
    this.max_x = 0;
    
    this.w = w;
    this.h = h;
    this.label_size = (this.h-40)/(this.tree_nodes*1.5);
    
    console.log(this.tree_len);
    console.log(this.tree_nodes);
    console.log(this.tree);
    this.svg = d3.select("#svg_div").append('svg')
      .attr('id', 'circlegraph_svg')
      .attr('width', this.w)
      .attr('height', this.h)

    this.range = [20, this.w-right_buffer];
    this.domain = [0, this.tree_len];
    this.xscale = d3.scaleLinear().domain(this.domain).range(this.range);
    this.yscale = d3.scaleLinear().domain([0, this.tree_nodes]).range([20, this.h-20]);

    this.draw_tree(this.tree, this.tree_nodes/2, 0)
  }

  draw_tree(tree, y_center, x) {
    //console.log(tree, y_center, x);
    // recursive tree drawing
    if (tree.children) {
      let node_counts = tree.children.map(t => t.node_count);
      let new_y = y_center - tree.node_count/2;
      let new_centers = [];
      for (let i=0; i<node_counts.length; i++) {
        let new_center = new_y+node_counts[i]/2;
        new_centers.push(new_center);
        let new_x = x+tree.children[i].length;
        this.svg.append('line')
          .attr('x1', this.xscale(x))
          .attr('x2', this.xscale(new_x))
          .attr('y1', this.yscale(new_center))
          .attr('y2', this.yscale(new_center))
          .attr('stroke', 'black');
        this.draw_tree(tree.children[i], new_center, new_x);
        new_y += node_counts[i];
      }
      this.svg.append('line')
        .attr('x1', this.xscale(x))
        .attr('x2', this.xscale(x))
        .attr('y1', this.yscale(new_centers[0]))
        .attr('y2', this.yscale(new_centers[new_centers.length-1]))
        .attr('stroke', 'black');
    } else {
      this.strain_locations[tree.name] = y_center;
      this.max_x = Math.max(this.max_x, this.xscale(x));
      //console.log(tree.name);
      this.svg.append('circle')
        .datum(tree.name.replace('.', '_'))
        .attr('class', 'tree_node')
        .attr('id', 'tree_node_'+tree.name.replace('.', '_'))
        .attr('cx', this.xscale(x)+5)
        .attr('cy', this.yscale(y_center))
        .attr('r', 5)
        .attr('fill', 'black')
        .attr('stroke', 'none')
        .attr('stroke-width', 3) // for highlighting focal
        .attr('opacity', 0.7)
        .on('mouseover', function(e, d) {
          d3.select(this).attr('opacity', 0.9).attr('fill', 'red');
          d3.select('#panpath_'+d)
            .attr('stroke-width', hover_stroke_width)
            .attr('opacity', 1)
        }).on('mouseout', function(e, d) {
          d3.select(this).attr('opacity', 0.7).attr('fill', 'black');
          d3.select('#panpath_'+d)
            .attr('stroke-width', 0.25)
            .attr('opacity', (core_only || show_paths_always) ? 0.5 : 0)
        }).on('click', function(e, d) {
          d3.selectAll('.tree_node').attr('stroke', 'none');
          remake_graph(d);
        });
        //console.log(this.label_size)
        this.svg.append('text')
          .style('cursor', 'default')
          .attr('x', this.xscale(x)+15)
          .attr('y', this.yscale(y_center)+this.label_size/3)
          .attr('font-size', this.label_size)
          .attr('fill', 'black')
          .text(tree.name);
    }
  }

  draw_junction_inset(j) {
    console.log(j);

    this.svg.selectAll('.junction_inset_subblock').remove();
    this.svg.selectAll('.junction_inset_subtext').remove();
    this.svg.select('#junction_inset_x_axis').remove();
    if (j) {
      const junc_xscale = d3.scaleLinear()
        .domain([0, j.maxLen])
        .range([this.max_x+100, this.w-20])

      const xAxis = d3.axisBottom(junc_xscale)
        .ticks(10, "~s");  // Adjust tick formatting as needed

      // Append x-axis to the SVG
      this.svg.append("g")
        .attr('id', 'junction_inset_x_axis')
        .attr("transform", `translate(0, ${this.h-20})`)  // Position the x-axis
        .call(xAxis);

      for (let strain of Object.keys(this.strain_locations)) {
        let top = this.yscale(this.strain_locations[strain]-0.4);
        let h = this.yscale(this.strain_locations[strain]+0.4)-this.yscale(this.strain_locations[strain]-0.4);
        if (j[strain].length == 1) { // structural change, no junction
          this.svg.append('text')
            .attr('class', 'junction_inset_subtext')
            .attr('x', junc_xscale(0))
            .attr('y', top+h)
            .attr('font-size', h*1.3)
            .attr('fill', 'black')
            .text(j[strain][0])
        } else {
          let running_sum = 0;
          for (let block_id of j[strain]) {
            let left = junc_xscale(running_sum);
            this.svg.append('rect')
              .datum(circlegraph.block_hash[block_id])
              .attr('class', 'junction_inset_subblock')
              .attr('x', left)
              .attr('width', (d) => junc_xscale(d.sequence.length+running_sum)-left)
              .attr('y', top)
              .attr('height', h)
              .attr('fill', circlegraph.colormap[block_id])
              .attr('stroke', '#FFF')
              .attr('stroke-width', 0.5)
            running_sum += circlegraph.block_hash[block_id].sequence.length;
          }
        }
      }
      d3.selectAll('.junction_inset_subblock')
        .on('mouseover', (e, d) => {
          d3.selectAll('.junction_indicator')
            .attr('opacity', (td) => d.in_junctions.indexOf(td.junction)>-1 ? 1 : 0);
        })
      if (annotated) {
        d3.selectAll('.junction_inset_subblock')
          .on('mousemove', (e, d) => {
            show_tooltip(e.offsetX, e.offsetY, d.annotations.join('\n'))
          })
          .on('mouseout', (e, d) => {
            d3.selectAll('.junction_indicator').attr('opacity', 0)
            hide_tooltip();
          })
      }
    }

  }
}

class CircleOccupancy {

  constructor(radial_chunks, chunk_buffer) {
    this.rad_chunks = radial_chunks;
    this.chunk_buf = chunk_buffer;
    this.rings = [];
  }

  place_arc(left_theta, right_theta) {
    // the arc runs clockwise from left_theta to right_theta
    // Adding this.rad_chunks here to avoid negatives, then taking the modulus 
    // in case we are over the number of chunks
    const left_pos = (this.rad_chunks + Math.floor((left_theta*this.rad_chunks)/(2*Math.PI))-this.chunk_buf) % this.rad_chunks;
    const right_pos = (this.rad_chunks + Math.ceil((right_theta*this.rad_chunks)/(2*Math.PI))+this.chunk_buf+1) % this.rad_chunks;
    let fit_ring;
    let new_ring;
    for (let ring_count=0; ring_count<=this.rings.length; ring_count++) {
      // copying the ring to fill it in, if it fits we'll keep the copy
      // If we've reached the end of the list, make a fresh ring
      let copied_occupancy = (ring_count == this.rings.length) ? new Array(this.rad_chunks).fill(false) : this.rings[ring_count].slice();
      let chunk = left_pos;
      let fit = true;
      while (true) {
        if (copied_occupancy[chunk]) {
          fit = false;
          break;
        } else {
          copied_occupancy[chunk] = true;
        }
        chunk += 1;
        if (chunk == this.rad_chunks) chunk = 0; // around the origin
        if (chunk == right_pos) break
      }
      if (fit) {
        new_ring = copied_occupancy;
        fit_ring = ring_count;
        break;
      }
    }
    if (fit_ring == this.rings.length) {
      this.rings.push(new_ring);
    } else {
      this.rings[fit_ring] = new_ring;
    }
    return fit_ring;
  }
}

// MAIN CLASS
class CircleGraph {

  constructor(pangraph_json_data) {
    console.log('constructing...');
    console.log(pangraph_json_data);

    this.focal_r = 160;
    this.block_h = 20;
    this.off_block_h = 8;
    this.c = 250;
    this.juncplot_left = 600
    this.juncplot_w = 250;
    this.juncplot_top = this.c-this.juncplot_w/2;
    this.juncplot_h = this.juncplot_w;
    this.junction_block_max_h = 80;

    this.data = pangraph_json_data;
    this.nStrains = this.data.paths.length;
    this.block_hash = this.data.blocks.reduce((map, block) => {
      block.is_core = this.block_is_core(block.id);
      block.not_dup = this.block_is_not_dup(block.id);
      block.in_junctions = [];
      map[block.id] = block;
      return map
    }, {});
    
    // DEPRECATED if (Object.keys(this.data.paths[0].blocks[0]).indexOf('annotated_seqs')>-1) {
    //  this.get_block_annotations();
    if (Object.keys(this.data.blocks[0]).indexOf('annotations')>-1) {
      annotated = true;
      console.log(this.block_hash);
    } else {
      annotated = false;
    }

    // assigning colors
    this.colormap = {};
    let counter = 0
    for (let p of this.data.paths) {
      for (let b of p.blocks) {
        this.colormap[b.id] = seaborn_colorblind[counter % seaborn_colorblind.length];
        counter += 1;
      }
    }

    this.block_positions = {};
    this.secondary_paths = {};
    //console.log(this.block_hash);
    this.w = 900;
    this.h = 500;
    this.svg = d3.select("#svg_div").append('svg')
      .attr('id', 'circlegraph_svg')
      .attr('width', this.w)
      .attr('height', this.h);
      //.style('background-color', 'pink');

    //d3.select("svg").call(d3.behavior.zoom());
    //this.setup_drag();
  }

  block_is_core(block_id) {
    const paths = this.data.paths.filter(p => 
      p.blocks.filter(b => b.id == block_id).length==1
    );
    return paths.length == this.data.paths.length;
  }

  block_is_not_dup(block_id) {
    const paths = this.data.paths.filter(p => 
      p.blocks.filter(b => b.id == block_id).length>1
    );
    return paths.length == 0;
  }

  /*
  DEPRECATED
  get_block_annotations() {
    for (let p of this.data.paths) {
      let name = p.name;
      for (let b of p.blocks) {
        let bid = b.id;
        let block_info = this.block_hash[bid];
        block_info.annotations = {};
        if (b.annotated_seqs) {
          for (let a of Object.values(b.annotated_seqs)) {
            // TODO should use locus id here as a key
            let gene_key = '<p><b>'+a.gene+'</b>: '+a.annotation+'</p>';
            if (Object.keys(block_info.annotations).indexOf(gene_key)==-1) {
              // new annotation
              block_info.annotations[gene_key] = a;
            }
          }
        }
      }
    }
  }
  */

  setup_drag() {
    const self = this;
    this.drag_r = this.focal_r-50;

    this.drag_el = this.svg.append('g');

    this.drag_circle = this.drag_el.append('circle')
      .attr('cx', self.c)
      .attr('cy', self.c)
      .attr('r', this.drag_r)
      .attr('stroke', '#333')
      .attr('fill', 'none')
      .attr('stroke-width', 30);

    this.drag_arc = this.drag_el.append('path')
      .attr('fill', 'none')
      .attr('stroke', 'red')
      .attr('opacity', 0.8)
      .attr('stroke-width', 30);

    // I realized D3 has an arc generator after I coded it up myself
    // for the alignments...
    const arcGenerator = d3.arc()
      .innerRadius(this.focal_r-50)
      .outerRadius(this.focal_r-50); 

    const dragBehavior = d3.drag()
      .on("start", dragStart)
      .on("drag", drag)
      .on("end", dragEnd);

    this.drag_el.call(dragBehavior);

    let startAngle;
    let clockwise;

    function getAngle(x, y) {
      const at2 = Math.atan2(y, x);
      if (at2 < 0) {
        return (Math.PI*2+at2)
      } else {
        return at2
      }
    }

    function dragStart(e) {
      const [x, y] = [e.x, e.y];
      startAngle = getAngle(x - self.c, y - self.c);
      clockwise = ['not sure yet', true];
    }

    function drag(e) {
      const [x, y] = [e.x, e.y];
      const currentAngle = getAngle(x - self.c, y - self.c);
      if (clockwise[0] === 'not sure yet') {
        // edge case
        let useAngle;
        if (Math.abs(currentAngle-startAngle)>Math.PI) {
          useAngle = (currentAngle > startAngle) ? currentAngle-Math.PI*2 : currentAngle+Math.PI*2;
        } else {
          useAngle = currentAngle;
        }
        if (useAngle > startAngle+Math.PI/10) { // arbitrary condition
          clockwise = ['locked_in', 1];
        } else if (useAngle < startAngle-Math.PI/10) {
          clockwise = ['locked_in', 0];
        } else {
          clockwise = ['not sure yet', useAngle > startAngle ? 1 : 0]
        }
      }
      let angle_difference;
      if ((clockwise[1] == 1) & (currentAngle < startAngle)) {
        angle_difference = Math.abs(Math.PI*2+currentAngle-startAngle);
      } else if ((clockwise[1] == 0) & (currentAngle > startAngle)) {
        angle_difference = Math.abs(Math.PI*-2+currentAngle-startAngle);
      } else {
        angle_difference = Math.abs(currentAngle-startAngle);
      }
      const large_arc = (angle_difference > Math.PI) ? 1 : 0;
      const arcPath = self.make_arc(startAngle, currentAngle, self.drag_r, clockwise[1], large_arc);
      // Update the arc path element (assuming you have one)
      self.drag_arc.attr("d", arcPath);
    }

    function dragEnd(e) {
      const [x, y] = [e.x, e.y];
      const currentAngle = getAngle(x - self.c, y - self.c);
      let start;
      let end;
      if ((clockwise[1] == 1) & (currentAngle < startAngle)) {
        start = startAngle;
        end = currentAngle;
      } else if ((clockwise[1] == 0) & (currentAngle > startAngle)) {
        start = currentAngle;
        end = startAngle;
      } else if (currentAngle > startAngle) {
        start = startAngle;
        end = currentAngle;
      } else {
        start = currentAngle;
        end = startAngle;
      }
      self.set_focal_region(Math.round(self.total_draw_len*start/(2*Math.PI)), 
                            Math.round(self.total_draw_len*end/(2*Math.PI)));
    }
  }

  get_nearby_core(path, position, direction) {
    // recursive function to find nearest core block
    // edge conditions
    if (direction == -1) {
      position = (position == -1) ? path.length-1 : position;
    } else {
      position = (position == path.length) ? 0 : position;
    }
    if (this.block_hash[path[position].id].is_core) {
      return path[position].id;
    } else {
      return this.get_nearby_core(path, position+direction, direction);
    }
  }

  process_path(path) {
    let all_redund_path_ids = new Set();
    const id_counts = path.reduce((id_counter, b) => {
        // counter one-liner using or shortcircuit
        id_counter[b.id] = (id_counter[b.id] || 0) + 1; 
        return id_counter;
      }, {});
    const repeated = Object.keys(id_counts).filter(id => id_counts[id]>1);
    let block_index = 0
    for (let b of path) {
      if (repeated.indexOf(b.id)>-1) {
        // label repeated blocks by the core blocks on either side
        let tmp_left = this.get_nearby_core(path, block_index, -1);
        let tmp_right = this.get_nearby_core(path, block_index, 1);
        // sorting the adjacent blocks to label the same regardless of direction
        let sorted_adjacents = [tmp_left, tmp_right].sort();
        let tmp_id = b.id+'_'+sorted_adjacents[0]+'_'+sorted_adjacents[0]+'_1';
        // if there are multiple repeats between core blocks,
        // number them
        let redundant_index = 2;
        while (all_redund_path_ids.has(tmp_id)) {
          tmp_id = tmp_id.split('_').concat([String(redundant_index)]).join('_');
          redundant_index += 1;
        }
        b.path_id = tmp_id;
        all_redund_path_ids.add(tmp_id);
      } else {
        b.path_id = b.id;
      }
      block_index += 1;
    }
  }

  get_block_in_and_out(block, block_info) {
    let strand_match = block_info.strand == block.strand;
    let direction = strand_match ? block_info.direction : -1*block_info.direction;
    let pos1 = (direction==1) ? block_info.l2 : block_info.l1;
    let pos2 = (direction==1) ? block_info.r2 : block_info.r1;
    let theta1 = (direction==1) ? block_info.left_theta : block_info.right_theta;
    let theta2 = (direction==1) ? block_info.right_theta : block_info.left_theta;
    return {
      'start': pos1, 
      'end': pos2, 
      'dir': direction, 
      'r': block_info.r,
      'h': block_info.h, 
      'theta1': theta1, 
      'theta2': theta2,
      'pid': block.path_id
    };
  }

  draw_path(raw_path, strain_id) {
    const path = draw_core_only ? raw_path.filter(b => this.block_hash[b.id].is_core) : raw_path;
    let svg_paths = this.svg.append('g')
      .attr('id', 'panpath_'+strain_id.replace('.', '_'))
      .attr('stroke-width', 0.25)
      .attr('opacity', (core_only || show_paths_always) ? 0.5 : 0)
    let current_direction = 0;
    let p = '';
    let position_info = []; //[pos1, pos2, direction, radius]
    for (let block of path) {
      let block_info = this.block_positions[block.path_id]
      position_info.push(this.get_block_in_and_out(block, block_info));
    }
    for (let b=0; b<path.length; b++) {
      let pi1 = position_info[b];
      let pi2 = position_info[(b+1)%path.length];
      if (pi1.dir != current_direction) {
        svg_paths.append('path')
          .attr('stroke', (current_direction==1) ? colors.forward_path: colors.rev_path)
          .attr('fill', 'none')
          .attr('d', p);
        p = '';
        current_direction = pi1.dir;
      }
      if (p=='') p += 'M '+String(pi1.start[0])+' '+String(pi1.start[1])+' ';
      // line across block
      let r = (pi1.dir == 1) ? pi1.r+pi1.h : pi1.r;
      let arc_tmp = 'A ' + String(r)+' '+String(r)+' 0 ';
      p += arc_tmp + '0 1 '+String(pi1.end[0])+' '+String(pi1.end[1])+' ';
      if (pi1.r == pi2.r) { // step nowhere or up
        if (pi1.theta2==pi2.theta1) {
          // same theta, no need for bezier
          p += 'L '+String(pi2.start[0])+' '+String(pi2.start[1])+' ';
        } else {
          //console.log(pi1.theta2-pi2.theta1)
          let extend1 = (pi1.dir == 1) ? 1+vert_bez : 1-vert_bez;
          let extend2 = (pi2.dir == 1) ? 1+vert_bez : 1-vert_bez;
          let bezier1 = extend_off_circle(this.c, this.c, pi1.end, extend1);
          let bezier2 = extend_off_circle(this.c, this.c, pi2.start, extend2);
          p += 'C '+String(bezier1[0])+' '+String(bezier1[1])+', '
          p += String(bezier2[0])+' '+String(bezier2[1])+', '
          p += String(pi2.start[0])+' '+String(pi2.start[1])+' '
        }
      } else if (pi1.r < pi2.r) {
        let bezier1 = extend_off_circle(this.c, this.c, pi1.end, 1+vert_bez);
        let bez_bump = 1.5*Math.abs(pi2.theta1-pi1.theta2)+0.001;
        let bezier2_theta = (pi2.dir == 1) ? pi2.theta1-bez_bump : pi2.theta1+bez_bump;
        let bezier2_r = (pi2.dir == 1) ? pi2.r+pi2.h : pi2.r;
        let bezier2 = point_on_circle(this.c, this.c, bezier2_r, bezier2_theta);
        p += 'C '+String(bezier1[0])+' '+String(bezier1[1])+', '
        p += String(bezier2[0])+' '+String(bezier2[1])+', '
        p += String(pi2.start[0])+' '+String(pi2.start[1])+' '
      } else { // step down
          let bezier2 = extend_off_circle(this.c, this.c, pi2.start, 1+vert_bez);
          let bez_bump = 1.5*Math.abs(pi2.theta1-pi1.theta2)+0.001;
          let bezier1_theta = (pi1.dir == 1) ? pi1.theta2+bez_bump : pi1.theta2-bez_bump;
          let bezier1_r = (pi1.dir == 1) ? pi1.r+pi1.h : pi1.r
          let bezier1 = point_on_circle(this.c, this.c, bezier1_r, bezier1_theta);
          p += 'C '+String(bezier1[0])+' '+String(bezier1[1])+', '
          p += String(bezier2[0])+' '+String(bezier2[1])+', '
          p += String(pi2.start[0])+' '+String(pi2.start[1])+' '
      }
      
    }
    svg_paths.append('path')
      .attr('stroke', (current_direction==1) ? colors.forward_path: colors.rev_path)
      .attr('fill', 'none')
      .attr('d', p);
  }

  get_anchor_theta(block_id, strand, start) {
    // strand is the strand of the secondary path block
    // which will be checked if it matches
    // start is true if we are looking for the start of the block
    // for the secondary path
    
    let block_info = this.block_positions[block_id];
    let strand_match = block_info.strand == strand;
    if (start) {
      return strand_match ? block_info.left_theta : block_info.right_theta;
    } else {
      return strand_match ? block_info.right_theta : block_info.left_theta;
    }
  }

  add_secondary_blocks(block_list) {
    const self = this;
    // the first and last blocks already exist, the others will be drawn
    let start = this.get_anchor_theta(block_list[0].path_id, block_list[0].strand, false);
    let end = this.get_anchor_theta(block_list[block_list.length-1].path_id, block_list[block_list.length-1].strand, true);
    let diff = end - start;
    let midpoint = (end + start)/2;
    let direction = this.block_positions[block_list[0].path_id].direction;
    //(diff > 0) ? 1 : -1;
    if (Math.abs(diff) > Math.PI) {
      midpoint = (midpoint + Math.PI) % (2*Math.PI);
      direction = (diff < 0) ? 1 : -1;
    }
    let insert_len = block_list.slice(1, block_list.length-1).reduce((summer, block) => {
      summer += this.block_hash[block.id.split('_')[0]].sequence.length;
      return summer;
    }, 0);
    let insert_angle_range = 2*Math.PI * insert_len / this.total_draw_len;
    let insert_left_theta = midpoint - insert_angle_range/2;
    let insert_right_theta = midpoint + insert_angle_range/2;
    let ring = this.circleOccupancy.place_arc(insert_left_theta, insert_right_theta);
    let base_r = (this.off_block_h*1+vert_bez)*ring+this.focal_r+7*this.block_h/4;

    let running_theta_len = 0;
    for (let block of block_list.slice(1, block_list.length-1)) {
      let theta_len = (this.block_hash[block.id].sequence.length*2*Math.PI)/this.total_draw_len;
      let left_theta = (direction == 1) ? insert_left_theta + running_theta_len : insert_right_theta - running_theta_len - theta_len;
      let right_theta = left_theta + theta_len;
      running_theta_len += theta_len;
      this.svg.append('path')
        .datum(block.path_id)
        .attr('stroke', '#FFF')
        .attr('stroke-width', 0.1)
        .attr('fill', ['green', 'orange', 'purple'][Math.floor(Math.random()*3)])
        .attr('d', this.draw_block_path(block, base_r, this.off_block_h, left_theta, right_theta, direction))
        .on('mouseover', function(e, d)  {
          console.log(self.block_positions[d]);
        });
    }
  }

  draw_secondary_strain(strain_index) {
    const strain_id = this.data.paths[strain_index].name
    if (core_only) {
      this.secondary_paths[strain_index] = this.data.paths[strain_index].blocks.filter(b => this.block_hash[b.id].is_core);
    } else if (remove_dups) {
      this.secondary_paths[strain_index] = this.data.paths[strain_index].blocks.filter(b => this.block_hash[b.id].not_dup);
    }else {
      this.secondary_paths[strain_index] = this.data.paths[strain_index].blocks;
    }
    
    this.process_path(this.secondary_paths[strain_index]);

    let junctions = [];
    let current_junc = [];
    let block_index = 0;
    while (block_index < this.secondary_paths[strain_index].length) {
      let b_id = this.secondary_paths[strain_index][block_index].path_id;
      //console.log(b_id, this.focal_strands[b_id], (b_id in this.focal_strands));
      if (b_id in this.block_positions) {
        if (current_junc.length > 0) {
          current_junc.push(this.secondary_paths[strain_index][block_index]);
          junctions.push(current_junc);
          current_junc = [this.secondary_paths[strain_index][block_index]];
        } else {
          current_junc = [this.secondary_paths[strain_index][block_index]];
        }
      } else if (current_junc.length > 0) {
        current_junc.push(this.secondary_paths[strain_index][block_index]);
      }
      block_index += 1;
    }
    // end condition TODO check this
    block_index = 0;
    while (block_index < this.secondary_paths[strain_index].length) {
      let b_id = this.secondary_paths[strain_index][block_index].path_id;
      if (b_id in this.block_positions) {
        current_junc.push(this.secondary_paths[strain_index][block_index]);
        junctions.push(current_junc);
        current_junc = [this.secondary_paths[strain_index][block_index]];
        break
      } else {
        current_junc.push(this.secondary_paths[strain_index][block_index]);
      }
      block_index += 1;
    }
    
    this.secondary_blocks = this.svg.append('g');
    for (let seg of junctions) {
      if (seg.length > 2) {
        this.add_secondary_blocks(seg);
      }
    }

    this.draw_path(this.secondary_paths[strain_index], strain_id);  
  }

  get_junc(blocks, block_1, block_2) {
    // TODO check this
    const block_list = blocks.map(b => b.id);
    const block_id1 = block_1.id;
    const block_id2 = block_2.id;
    const start = block_list.indexOf(block_id1);
    const end = block_list.indexOf(block_id2);
    // checking for junction orientation
    if ((block_1.strand==block_list[start].strand) != (block_2.strand==block_list[end].strand)) {
      return ['inversion between blocks'];
    }
    const diff = end-start;
    let block_subset;
    if (Math.abs(diff) > block_list.length/2) {
      // goes around list origin
      if (diff > 0) { // counterclockwise
        block_subset = block_list.slice(0, start+1).reverse().concat(block_list.slice(end, block_list.length).reverse());
      } else {
        block_subset = block_list.slice(start, block_list.length).concat(block_list.slice(0, end+1));
      }
    } else {
      if (diff > 0) {
        block_subset = block_list.slice(start, (block_list.length+end+1) % block_list.length);
      } else {
        block_subset = block_list.slice(end, (block_list.length+start+1) % block_list.length).reverse();
      }
    }
    if (block_subset.filter(b => this.block_hash[b].is_core).length > 2) {
      return ['structural change between blocks']
    }
    
    return block_subset;
  }

  get_junctions() {
    // TODO DOES NOT WORRY ABOUT STRAND RN
    const fpath = this.focal_path.filter(b => this.block_hash[b.id].is_core);
    this.junctions = {};
    for (let i=0; i<fpath.length; i++) {
      let j = (i+1) % fpath.length;
      this.junctions[fpath[i].id+'_'+fpath[j].id] = {'junction': fpath[i].id+'_'+fpath[j].id};
      let all_juncs = []
      let accessory_blocks = new Set()
      let maxLen = 0
      for (let p of this.data.paths) {
        let junc = this.get_junc(p.blocks, fpath[i], fpath[j]);
        this.junctions[fpath[i].id+'_'+fpath[j].id][p.name] = junc;
        all_juncs.push(junc.join(';'));
        if (junc.length > 2) {
          for (let bid of junc.slice(1,junc.length-1)) {
            accessory_blocks.add(bid);
            this.block_hash[bid].in_junctions.push(fpath[i].id+'_'+fpath[j].id);
          }
          let junc_len = junc.reduce((summer, b) => {
            // using consensus seq lengths
            summer += this.block_hash[b].sequence.length;
            return summer;
          }, 0);
          maxLen = Math.max(maxLen, junc_len);
        }
      }
      this.junctions[fpath[i].id+'_'+fpath[j].id]['nUniqueJunctions'] = new Set(all_juncs).size;
      this.junctions[fpath[i].id+'_'+fpath[j].id]['pangenomeSize'] = Array(...accessory_blocks).reduce((summer, b) => {
        // using consensus seq lengths
        summer += this.block_hash[b].sequence.length;
        return summer;
      }, 0);
      this.junctions[fpath[i].id+'_'+fpath[j].id]['maxLen'] = maxLen;
    }
    console.log(this.junctions);
    this.max_nUniqueJunctions = Math.max(...Object.values(this.junctions).map(j => j.nUniqueJunctions));
    this.max_pangenomeSize = Math.max(...Object.values(this.junctions).map(j => j.pangenomeSize));
    console.log(this.max_nUniqueJunctions, this.max_pangenomeSize);
    this.draw_junction_scatterplot();
    this.draw_junction_stat_blocks();
  }

  draw_junction_scatterplot() {
    const self = this;
    this.juncplot_xscale = d3.scaleLog()
      .domain([0.9, this.max_pangenomeSize*2])
      .range([this.juncplot_left, this.juncplot_left+this.juncplot_w])
    this.juncplot_yscale = d3.scaleLog()
      .domain([0.9, this.max_nUniqueJunctions*1.05])
      .range([this.juncplot_top+this.juncplot_h, this.juncplot_top])

    // Create the x-axis
    const xAxis = d3.axisBottom(this.juncplot_xscale)
        .ticks(10, "~s");  // Adjust tick formatting as needed

    // Create the y-axis
    const yAxis = d3.axisLeft(this.juncplot_yscale)
        .ticks(10, "~s");  // Adjust tick formatting as needed

    // Append x-axis to the SVG
    this.svg.append("g")
        .attr("transform", `translate(0, ${this.juncplot_top + this.juncplot_h})`)  // Position the x-axis
        .call(xAxis);

    // Append y-axis to the SVG
    this.svg.append("g")
        .attr("transform", `translate(${this.juncplot_left}, 0)`)  // Position the y-axis
        .call(yAxis);

    // Add x-axis label
    this.svg.append("text")
        .attr("class", "x label")
        .attr("text-anchor", "middle")
        .attr("x", this.juncplot_left+this.juncplot_w/2)
        .attr("y", this.juncplot_top+this.juncplot_h + 40)  // Position below the x-axis
        .text("Local pangenome length (bp)");

    // Add y-axis label
    this.svg.append("text")
        .attr("class", "y label")
        .attr("text-anchor", "middle")
        .attr("x", this.juncplot_left)
        .attr("y", this.juncplot_top - 10)  
        .text("Number of unique junctions");  

    this.svg.selectAll('.juncplot_point')
      .data(Object.values(this.junctions))
      .enter()
      .append('circle')
        .attr('class', 'juncplot_point')
        .attr('cx', (d) => this.juncplot_xscale(d.pangenomeSize))
        .attr('cy', (d) => this.juncplot_yscale(d.nUniqueJunctions)+Math.random()*10) // JITTER
        .attr('r', 3)
        .attr('fill', '#333')
        .attr('stroke', 'red')
        .attr('stroke-width', 0)
        .on('mouseover', function(e, d)  {
          d3.select(this).attr('opacity', 0.8);
          d3.select(this).attr('stroke-width', 2);
          d3.selectAll('.junction_block').attr('stroke-width', (td) => (td.junction == d.junction) ? 4 : (td.junction==clickedJunction) ? 6 : 0.5);
          treeObj.draw_junction_inset(d);
        })
        .on('mouseout', function(e, d) {
          d3.select(this).attr('opacity', 1)
          d3.select(this).attr('stroke-width', (td) => (td.junction==clickedJunction) ? 4 : 0);
          d3.selectAll('.junction_block')
            .attr('stroke-width', (td) => (td.junction==clickedJunction) ? 6 : 0.5);
          treeObj.draw_junction_inset(self.junctions[clickedJunction]);
        })
        .on('click', function(e, d) {
          if (clickedJunction == d.junction) {
            clickedJunction = null;
          } else {
            clickedJunction = d.junction;
          }
          treeObj.draw_junction_inset(self.junctions[clickedJunction]);
          d3.selectAll('.juncplot_point').attr('stroke-width', (td) => (td.junction==clickedJunction) ? 4 : 0);
          d3.selectAll('.junction_block')
            .attr('stroke-width', (td) => (td.junction==clickedJunction) ? 6 : 0.5)
        });
  }

  draw_junction_stat_blocks() {
    // draws blocks onto the outside of the circle
    // representing junction stats - n unique junctions
    // and total pangenome length
    const self = this;
    for (let j of Object.values(this.junctions)) {
      let anchor_block_ids = j.junction.split('_')
      // TODO this is hacky, just reusing code and making sure teh strands match
      let start = this.get_anchor_theta(anchor_block_ids[0], this.block_positions[anchor_block_ids[0]].strand, false);
      let end = this.get_anchor_theta(anchor_block_ids[1], this.block_positions[anchor_block_ids[1]].strand, true);
      let diff = end - start;
      let midpoint = (end + start)/2;
      if (Math.abs(diff) > Math.PI) {
        midpoint = (midpoint + Math.PI) % (2*Math.PI);
      }
      let insert_len = j.pangenomeSize / this.nStrains;
      let h = this.junction_block_max_h * (j.nUniqueJunctions/this.max_nUniqueJunctions);
      let insert_angle_range = 2*Math.PI * insert_len / this.total_draw_len;
      let insert_left_theta = midpoint - insert_angle_range/2;
      let insert_right_theta = midpoint + insert_angle_range/2;

      this.svg.append('circle')
        .datum(j)
        .attr('class', 'junction_indicator')
        .attr('fill', 'black')
        .attr('r', 4)
        .attr('opacity', 0)
        .attr('cx', point_on_circle(this.c, this.c, this.focal_r+this.block_h+5, midpoint)[0])
        .attr('cy', point_on_circle(this.c, this.c, this.focal_r+this.block_h+5, midpoint)[1])

      this.svg.append('path')
        .datum(j)
        .attr('class', 'junction_block')
        .attr('stroke', 'red')
        .attr('stroke-width', 0.5)
        .attr('fill', ['green', 'orange', 'purple'][Math.floor(Math.random()*3)])
        .attr('d', this.draw_block_path(null, this.focal_r+this.block_h*1.5, h, insert_left_theta, insert_right_theta, null))
        .on('mouseover', function(e, d)  {
          console.log(d);
          d3.select(this).attr('opacity', 0.8);
          d3.select(this).attr('stroke-width', 4);
          d3.selectAll('.juncplot_point').attr('stroke-width', (td) => (td.junction == d.junction) ? 2 : (td.junction==clickedJunction) ? 4 : 0);
          treeObj.draw_junction_inset(d);
        })
        .on('mouseout', function(e, d) {
          d3.select(this).attr('opacity', 1)
          d3.select(this).attr('stroke-width', (td) => (td.junction==clickedJunction) ? 6 : 0.5);
          d3.selectAll('.juncplot_point').attr('stroke-width', (td) => (td.junction==clickedJunction) ? 4 : 0);
          treeObj.draw_junction_inset(self.junctions[clickedJunction]);
        })
        .on('click', function(e, d) {
          if (clickedJunction == d.junction) {
            clickedJunction = null;
          } else {
            clickedJunction = d.junction;
          }
          treeObj.draw_junction_inset(self.junctions[clickedJunction]);
          d3.selectAll('.junction_block').attr('stroke-width', (td) => (td.junction==clickedJunction) ? 6 : 0.5);
          d3.selectAll('.juncplot_point').attr('stroke-width', (td) => (td.junction==clickedJunction) ? 4 : 0)
        });;
    }
  }

  draw_focal_strain(strain_index) {
    const focal_strain_id = this.data.paths[strain_index].name;
    d3.select('.tree_node').attr('stroke', 'none');
    d3.select('#tree_node_'+focal_strain_id.replace('.', '_')).attr('stroke', '#F80');
    // setting up the circleOccupancy object
    this.circleOccupancy = new CircleOccupancy(2000, 1);
    // copying the path out and added sequence lengths
    // TODO I shouldn't need to store these paths twice
    this.focal_path = this.data.paths[strain_index].blocks;
    if (core_only) {
      this.focal_path_use = this.focal_path.filter(b => this.block_hash[b.id].is_core);
    } else if (remove_dups) {
      this.focal_path_use = this.focal_path.filter(b => this.block_hash[b.id].not_dup);
    } else {
      this.focal_path_use = this.focal_path;
    }
    
    this.process_path(this.focal_path_use);
    this.total_block_len = this.focal_path_use.reduce((summer, block) => summer + this.block_hash[block.id].sequence.length, 0);
    this.total_draw_len = this.total_block_len / (1-path_view_frac);
    this.single_path_len = (this.total_draw_len * path_view_frac) / this.focal_path_use.length;
    // dividing the fraction of circumf allowed for paths / the # of blocks
    this.block_spacing = (2 * Math.PI * path_view_frac) / this.focal_path_use.length;
    let running_start = 0;
    let block_index = 0
    this.focal_strands = {};
    for (let b of this.focal_path_use) {
      b.start = running_start;
      b.seq_len = this.block_hash[b.id].sequence.length
      b.end = b.start + b.seq_len;
      running_start += b.seq_len + this.single_path_len;
      //b.svg_position = this.get_anchors(b);
      this.focal_strands[b.path_id] = b.strand;
      this.draw_focal_block(b, greys[block_index % greys.length])
      block_index += 1
    }
    console.log(this.focal_path_use);
    //this.draw_path(this.focal_path);
    this.get_junctions();
  }

  draw_focal_block(b, fill_c) {
    this.svg.append('path')
      .attr('class', 'focal_block')
      .attr('stroke', '#FFF')
      .attr('stroke-width', 0.05)
      .attr('fill', fill_c)
      .attr('d', this.make_focal_path(b));
  }

  get_anchors(d) {
    const start_theta = 2*Math.PI*(d.start/this.total_draw_len);
    const end_theta = 2*Math.PI*(d.end/this.total_draw_len);
    //const top = point_on_circle(this.c, this.c, this.focal_r+this.block_h/2, mid_theta);
    //const bottom = point_on_circle(this.c, this.c, this.focal_r-this.block_h/2, mid_theta);
    return [start_theta, end_theta];
  }

  make_focal_path(block) {
    const left_theta = 2*Math.PI*(block.start/this.total_draw_len);
    // modulo makes sure we are in [0, 2pi)
    const right_theta = 2*Math.PI*(block.end/this.total_draw_len) % (Math.PI*2);
    return this.draw_block_path(block, this.focal_r, this.block_h, left_theta, right_theta, 1);
  }

  draw_block_path(block, r, h, left_theta, right_theta, direction) {
    let p = ''
    const arc_outside = 'A ' + String(r+h)+' '+String(r+h)+' 0 ';
    const arc_inside = 'A ' + String(r)+' '+String(r)+' 0 ';
    let s1 = point_on_circle(this.c, this.c, r, left_theta);
    let s2 = point_on_circle(this.c, this.c, r+h, left_theta);
    let e1 = point_on_circle(this.c, this.c, r, right_theta);
    let e2 = point_on_circle(this.c, this.c, r+h, right_theta);
    p += 'M'+String(s1[0])+' '+String(s1[1])+' '
    p += 'L'+String(s2[0])+' '+String(s2[1])+' '
    p += arc_outside + '0 1 '+String(e2[0])+' '+String(e2[1])+' '
    p += 'L'+String(e1[0])+' '+String(e1[1])+' '
    p += arc_inside + '0 0 '+String(s1[0])+' '+String(s1[1])+''
    if (block) {
      this.block_positions[block.path_id] = {
        'path_id': block.path_id,
        'strand': block.strand,
        'direction': direction,
        'left_theta': left_theta,
        'right_theta': right_theta,
        'r': r,
        'h': h,
        'l1': s1, 'l2': s2, 'r1': e1, 'r2': e2
      }
    }
    return p;
  }
}

function remake_graph(strain_name) {
  let strain_index = circlegraph.data.paths.map(p => p.name.replace('.', '_')).indexOf(strain_name);
  circlegraph.svg.remove();
  circlegraph = new CircleGraph(graph_data);
  circlegraph.draw_focal_strain(strain_index);
  //cc.draw_secondary_strain(9);
  for (let i=0; i<circlegraph.data.paths.length; i++) {
    circlegraph.draw_secondary_strain(i);
  }
}

function make_graph(strain_index) {
  circlegraph = new CircleGraph(graph_data);
  circlegraph.draw_focal_strain(strain_index);
  //cc.draw_secondary_strain(9);
  for (let i=0; i<circlegraph.data.paths.length; i++) {
    circlegraph.draw_secondary_strain(i);
  }
  setup_tooltip();
}

function load_example(ex) {
  console.log(ex);
  d3.selectAll('.input_stuff').remove();
  d3.select('#form_div').append('h3')
    .html('Example data: '+ex);
  d3.text(dsets[ex].treefile).then(function(tree_data) {
    treeObj = new SimpleTree(tree_data, 900, 400, 500);
    d3.json(dsets[ex].graphfile).then(function(data) {
      // Log the data to the console to verify it is loaded correctly
      graph_data = data;
      make_graph(graph_data.paths.length-1);
    }).catch(function(error) {
        console.log(error);
    });
  }).catch(function(error) {
    console.log(error);
  });
}

function goforit() {
  if (use_test_data) {
    load_example(use);
  } else {
    const graph_file_div = d3.select('#form_div').append('div');
    graph_file_div.append('label')
      .attr('class', 'file_in_label input_stuff')
      .attr('for', 'graph_file_input')
      .html('Load pangraph json file: ')
    graph_file_div.append('input')
      .attr('class', 'file_in input_stuff')
      .attr('type', 'file')
      .attr('id', 'graph_file_input')
      .on('input', load_graph_file)

    let button_holder = d3.select('#form_div').append('div')
      .attr('id', 'button_holder')
      .attr('class', 'input_stuff');

    button_holder.append('h3')
      .html('Or look at an example dataset:');
      
    button_holder.selectAll('.example_button')
      .data(dset_names)
      .enter()
        .append('div')
          .attr('class', 'example_button')
          .html((d) => d)
          .on('click', (e, d) => {
            load_example(d)
          });

  }
}

function load_graph_file() {
  d3.select('#button_holder').remove();
  let file = document.getElementById('graph_file_input').files[0]
  let reader = new FileReader();
    reader.readAsText(file);
    reader.onload = function() {
      graph_data = JSON.parse(reader.result);
      console.log(graph_data);
      const tree_file_div = d3.select('#form_div').append('div');
      tree_file_div.append('label')
        .attr('class', 'file_in_label input_stuff')
        .attr('for', 'tree_file_input')
        .html('Load core tree newick file: ')
      tree_file_div.append('input')
        .attr('class', 'file_in input_stuff')
        .attr('type', 'file')
        .attr('id', 'tree_file_input')
        .on('input', load_tree_file);
      d3.select('#form_div').append('div')
        .attr('class', 'example_button input_stuff')
        .attr('id', 'load_wo_tree_button')
        .html('Load without a tree')
        .on('click', load_without_tree);
    }
}

function load_without_tree() {
  d3.select('#load_wo_tree_button').remove();
  const tree_data = '(' + graph_data.paths.map((p) => p.name+':1').join(',') + ')';
  treeObj = new SimpleTree(tree_data, 900, 400, 800);
  make_graph(graph_data.paths.length-1);
}

function load_tree_file() {
  d3.select('#load_wo_tree_button').remove();
  let file = document.getElementById('tree_file_input').files[0]
  let reader = new FileReader();
    reader.readAsText(file);
    reader.onload = function() {
      const tree_data = reader.result;
      treeObj = new SimpleTree(tree_data, 900, 400, 500);
      make_graph(graph_data.paths.length-1);
    }
}

function start_from_streamlit() {
  treeObj = new SimpleTree(pangraphVizData.tree_string, 900, 400, 500);
  graph_data = pangraphVizData.graph_json;
  make_graph(graph_data.paths.length-1);
}