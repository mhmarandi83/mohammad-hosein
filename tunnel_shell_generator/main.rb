# tunnel_shell_generator/main.rb
require 'sketchup.rb'
require 'json'

module TunnelShellGenerator

  module Geometry
    def self.interpolate_points(p1, p2, t)
      Geom::Point3d.linear_combination(1.0 - t, p1, t, p2)
    end

    def self.project_point_to_segment(p, a, b)
      ab = b - a
      ab_len = ab.length
      return a if ab_len < 0.00001
      ab_unit = ab.normalize
      ap = p - a
      projection = ap.dot(ab_unit)
      if projection <= 0
        a
      elsif projection >= ab_len
        b
      else
        a.offset(ab_unit, projection)
      end
    end
  end

  module CurveMatching
    def self.extract_curves(selection)
      edges = selection.grep(Sketchup::Edge)
      return nil if edges.empty?
      chains = group_connected_edges(edges)
      return nil if chains.length != 2
      chains
    end

    def self.group_connected_edges(edges)
      unvisited = edges.dup
      chains = []
      until unvisited.empty?
        chain = [unvisited.shift]
        loop do
          grown = false
          endpoints = get_endpoints(chain)
          v_start = endpoints.first
          v_end = endpoints.last
          break if v_start.nil? || v_end.nil?
          unvisited.each_with_index do |edge, idx|
            if edge.start == v_end || edge.end == v_end
              chain.push(edge)
              unvisited.delete_at(idx)
              grown = true
              break
            elsif edge.start == v_start || edge.end == v_start
              chain.unshift(edge)
              unvisited.delete_at(idx)
              grown = true
              break
            end
          end
          break unless grown
        end
        chains << sort_edge_chain(chain)
      end
      chains
    end

    def self.get_endpoints(edges)
      vertex_counts = Hash.new(0)
      edges.each do |e|
        vertex_counts[e.start] += 1
        vertex_counts[e.end] += 1
      end
      endpoints = vertex_counts.select { |_, count| count == 1 }.keys
      if endpoints.length == 2
        endpoints
      else
        [edges.first.start, edges.last.end]
      end
    end

    def self.sort_edge_chain(edges)
      return edges if edges.length <= 1
      endpoints = get_endpoints(edges)
      current_vertex = endpoints.first || edges.first.start
      sorted = []
      remaining = edges.dup
      until remaining.empty?
        next_edge = remaining.find { |e| e.start == current_vertex || e.end == current_vertex }
        break unless next_edge
        sorted << next_edge
        remaining.delete(next_edge)
        current_vertex = (next_edge.start == current_vertex) ? next_edge.end : next_edge.start
      end
      sorted
    end

    class Polyline
      attr_reader :points, :lengths, :total_length

      def initialize(vertices)
        @points = vertices.map { |v| v.respond_to?(:position) ? v.position : v }
        @points.uniq!
        calculate_lengths
      end

      def calculate_lengths
        @lengths = [0.0]
        accumulated = 0.0
        (0...@points.length - 1).each do |i|
          dist = @points[i].distance(@points[i + 1])
          accumulated += dist
          @lengths << accumulated
        end
        @total_length = accumulated
      end

      def get_point_at(t)
        return @points.first if t <= 0.0
        return @points.last if t >= 1.0
        target_len = t * @total_length
        idx = nil
        @lengths.each_with_index do |x, i|
          if x >= target_len
            idx = i
            break
          end
        end
        return @points.last if idx.nil?
        return @points[idx] if idx == 0
        l_prev = @lengths[idx - 1]
        l_next = @lengths[idx]
        seg_t = (target_len - l_prev) / (l_next - l_prev)
        Geometry.interpolate_points(@points[idx - 1], @points[idx], seg_t)
      end

      def project_point(p)
        min_dist = Float::INFINITY
        closest_pt = nil
        (0...@points.length - 1).each do |i|
          a = @points[i]
          b = @points[i + 1]
          proj = Geometry.project_point_to_segment(p, a, b)
          d = p.distance(proj)
          if d < min_dist
            min_dist = d
            closest_pt = proj
          end
        end
        closest_pt || p
      end
    end

    def self.build_polyline(edges)
      vertices = []
      if edges.length == 1
        vertices = [edges.first.start, edges.first.end]
      else
        endpoints = get_endpoints(edges)
        current_vertex = endpoints.first || edges.first.start
        vertices << current_vertex
        edges.each do |e|
          next_vertex = (e.start == current_vertex) ? e.end : e.start
          vertices << next_vertex
          current_vertex = next_vertex
        end
      end
      Polyline.new(vertices)
    end

    def self.align_polylines(poly1, poly2)
      p1_start = poly1.points.first
      p1_end = poly1.points.last
      p2_start = poly2.points.first
      p2_end = poly2.points.last
      d_same = p1_start.distance(p2_start) + p1_end.distance(p2_end)
      d_opposite = p1_start.distance(p2_end) + p1_end.distance(p2_start)
      if d_opposite < d_same
        reversed_points = poly2.points.reverse
        poly2 = Polyline.new(reversed_points)
      end
      [poly1, poly2]
    end

    def self.match_points(poly1, poly2, samples_t)
      left_pts = []
      right_pts = []
      samples_t.each do |t|
        p1 = poly1.get_point_at(t)
        p2 = poly2.project_point(p1)
        left_pts << p1
        right_pts << p2
      end
      { left: left_pts, right: right_pts }
    end
  end

  module Sampling
    def self.sample_curve(polyline, params)
      step_dist = params[:sampling_distance]
      adaptive = params[:adaptive_sampling]
      return uniform_sampling(polyline, step_dist) unless adaptive
      adaptive_sampling(polyline, step_dist)
    end

    def self.uniform_sampling(polyline, step_dist)
      total_len = polyline.total_length
      return [0.0, 1.0] if total_len <= step_dist
      num_steps = (total_len / step_dist).ceil
      samples = []
      (0..num_steps).each do |i|
        samples << i.to_f / num_steps
      end
      samples
    end

    def self.adaptive_sampling(polyline, base_step)
      total_len = polyline.total_length
      return [0.0, 1.0] if total_len <= base_step
      t_samples = [0.0]
      current_len = 0.0
      min_step = base_step * 0.15
      max_step = base_step
      angle_threshold = 0.07
      current_step = base_step
      while current_len < total_len
        next_len = current_len + current_step
        if next_len >= total_len
          t_samples << 1.0
          break
        end
        t_curr = current_len / total_len
        t_next = next_len / total_len
        t_mid = (t_curr + t_next) / 2.0
        p_curr = polyline.get_point_at(t_curr)
        p_mid = polyline.get_point_at(t_mid)
        p_next = polyline.get_point_at(t_next)
        v1 = (p_mid - p_curr).normalize
        v2 = (p_next - p_mid).normalize
        dot = v1.dot(v2)
        dot = [-1.0, [1.0, dot].min].max
        angle = Math.acos(dot)
        if angle > angle_threshold && current_step > min_step
          current_step = [current_step * 0.5, min_step].max
        else
          current_len = next_len
          t_samples << (current_len / total_len)
          if angle < (angle_threshold * 0.3)
            current_step = [current_step * 1.5, max_step].max
          end
        end
      end
      t_samples.uniq.sort
    end
  end

  module FrameGenerator
    class Frame
      attr_reader :center, :tangent, :up, :right, :width
      def initialize(center, tangent, up, right, width)
        @center = center
        @tangent = tangent
        @up = up
        @right = right
        @width = width
      end
    end

    def self.generate_frames(left_points, right_points)
      n = left_points.length
      frames = []
      centers = []
      (0...n).each do |i|
        centers << Geom::Point3d.new(
          (left_points[i].x + right_points[i].x) / 2.0,
          (left_points[i].y + right_points[i].y) / 2.0,
          (left_points[i].z + right_points[i].z) / 2.0
        )
      end
      
      tangents = []
      (0...n).each do |i|
        if i == 0
          t = centers[1] - centers[0]
        elsif i == n - 1
          t = centers[n - 1] - centers[n - 2]
        else
          t = centers[i + 1] - centers[i - 1]
        end
        t = t.normalize
        tangents << t
      end
      
      t0 = tangents[0]
      world_up = Geom::Vector3d.new(0, 0, 1)
      if t0.parallel?(world_up)
        world_up = Geom::Vector3d.new(1, 0, 0)
      end
      
      initial_right = (right_points[0] - left_points[0]).normalize
      initial_up = (t0 * initial_right).normalize
      initial_right = (initial_up * t0).normalize
      
      up_vectors = [initial_up]
      right_vectors = [initial_right]
      
      (0...n - 1).each do |i|
        x_i = centers[i]
        x_next = centers[i + 1]
        t_i = tangents[i]
        t_next = tangents[i + 1]
        u_i = up_vectors[i]
        
        v1 = x_next - x_i
        c1 = v1.dot(v1)
        
        if c1 > 0.00001
          u_i_l = u_i - v1 * (2.0 * (v1.dot(u_i)) / c1)
          t_i_l = t_i - v1 * (2.0 * (v1.dot(t_i)) / c1)
          
          v2 = t_next - t_i_l
          c2 = v2.dot(v2)
          
          u_next = if c2 > 0.00001
                     u_i_l - v2 * (2.0 * (v2.dot(u_i_l)) / c2)
                   else
                     u_i_l
                   end
        else
          u_next = u_i
        end
        
        u_next = u_next.normalize
        r_next = (u_next * t_next).normalize
        
        up_vectors << u_next
        right_vectors << r_next
      end
      
      (0...n).each do |i|
        exact_right = (right_points[i] - left_points[i]).normalize
        exact_up = (tangents[i] * exact_right).normalize
        exact_right = (exact_up * tangents[i]).normalize
        width = left_points[i].distance(right_points[i])
        frames << Frame.new(centers[i], tangents[i], exact_up, exact_right, width)
      end
      frames
    end
  end

  module SectionGenerator
    class Section
      attr_reader :points, :width, :frame, :left_base, :right_base, :left_wall_top, :right_wall_top
      def initialize(points, width, frame, left_base, right_base, left_wall_top, right_wall_top)
        @points = points
        @width = width
        @frame = frame
        @left_base = left_base
        @right_base = right_base
        @left_wall_top = left_wall_top
        @right_wall_top = right_wall_top
      end
    end

    def self.generate_sections(frames, params)
      sections = []
      wall_height = params[:wall_height]
      roof_radius = params[:roof_radius]
      roof_segments = params[:roof_segments]
      gen_roof = params[:generate_roof]
      gen_walls = params[:generate_walls]
      gen_floor = params[:generate_floor]

      frames.each_with_index do |frame, idx|
        w = frame.width
        c = frame.center
        r = frame.right
        u = frame.up
        half_w = w / 2.0
        
        if gen_roof && roof_radius < half_w
          raise StandardError, "Radius is too small! At section #{idx + 1}, width is #{sprintf('%.2f', w.to_m)}m, which requires a minimum roof radius of #{sprintf('%.2f', half_w.to_m)}m. Requested radius was #{sprintf('%.2f', roof_radius.to_m)}m."
        end
        
        pts = []
        left_base = c.offset(r, -half_w)
        right_base = c.offset(r, half_w)
        left_wall_top = left_base.offset(u, wall_height)
        right_wall_top = right_base.offset(u, wall_height)

        if gen_floor
          pts << left_base
        end
        
        if gen_walls
          pts << left_base unless gen_floor
          pts << left_wall_top
        else
          pts << left_base if !gen_floor && !gen_roof
        end
        
        if gen_roof
          m = c.offset(u, wall_height)
          d = Math.sqrt(roof_radius**2 - half_w**2)
          arc_center = m.offset(u, -d)
          half_angle = Math.asin(half_w / roof_radius)
          (0..roof_segments).each do |s|
            ratio = s.to_f / roof_segments
            angle = -half_angle + (2.0 * half_angle * ratio)
            arc_r_offset = r.x * (roof_radius * Math.sin(angle))
            arc_u_offset = u.x * (roof_radius * Math.cos(angle))
            arc_pt = arc_center.offset(r, (roof_radius * Math.sin(angle)))
            arc_pt = arc_pt.offset(u, (roof_radius * Math.cos(angle)))
            if s == 0 && gen_walls
              next
            elsif s == roof_segments && gen_walls
              next
            end
            pts << arc_pt
          end
        end
        
        if gen_walls
          pts << right_wall_top
          pts << right_base
        else
          pts << right_base if !gen_floor && !gen_roof
        end
        
        if gen_floor
          pts << right_base unless pts.include?(right_base)
        end
        
        sections << Section.new(pts, w, frame, left_base, right_base, left_wall_top, right_wall_top)
      end
      sections
    end
  end

  module MeshBuilder
    def self.build(model, sections, params)
      num_sections = sections.length
      return if num_sections < 2
      pts_per_section = sections.first.points.length
      num_vertices = num_sections * pts_per_section
      num_polygons = (num_sections - 1) * (pts_per_section - 1)
      if params[:generate_end_caps]
        num_polygons += 2
      end
      mesh = Geom::PolygonMesh.new(num_vertices, num_polygons)
      vertex_indices = Array.new(num_sections) { Array.new(pts_per_section) }
      (0...num_sections).each do |i|
        (0...pts_per_section).each do |j|
          vertex_indices[i][j] = mesh.add_point(sections[i].points[j])
        end
      end
      (0...num_sections - 1).each do |i|
        (0...pts_per_section - 1).each do |j|
          v0 = vertex_indices[i][j]
          v1 = vertex_indices[i + 1][j]
          v2 = vertex_indices[i + 1][j + 1]
          v3 = vertex_indices[i][j + 1]
          mesh.add_polygon(v0, v1, v2, v3)
        end
      end
      if params[:generate_end_caps]
        start_indices = (0...pts_per_section).map { |j| vertex_indices[0][j] }.reverse
        mesh.add_polygon(*start_indices)
        end_indices = (0...pts_per_section).map { |j| vertex_indices[num_sections - 1][j] }
        mesh.add_polygon(*end_indices)
      end
      active_ents = model.active_entities
      group = active_ents.add_group
      group.name = "Tunnel Shell"
      tunnel_tag = model.layers.add("Tunnel Shell")
      group.layer = tunnel_tag
      group.entities.add_faces_from_mesh(mesh, true, nil)
      if params[:smooth_edges]
        smooth_and_soften_edges(group)
      end
      group
    end

    def self.smooth_and_soften_edges(group)
      threshold_angle = 35.0 * Math::PI / 180.0
      group.entities.grep(Sketchup::Edge).each do |edge|
        next if edge.faces.length != 2
        f1 = edge.faces[0]
        f2 = edge.faces[1]
        angle = f1.normal.angle_between(f2.normal)
        if angle < threshold_angle
          edge.soft = true
          edge.smooth = true
        end
      end
    end
  end

  module TunnelUI
    def self.show_dialog(callback)
      begin
        puts ">> Opening dialog..."
        dialog = UI::WebDialog.new(
          "Tunnel Shell Generator",
          true,
          true,
          350,
          580,
          200,
          200,
          true
        )
        dialog.set_html(get_html_content)

        dialog.add_action_callback("generate") do |_dialog, params_json|
          puts ">> Generate button clicked with params: #{params_json}"

          begin
            params_hash = JSON.parse(params_json)
            
            ruby_params = {
              roof_radius: params_hash['roof_radius'].to_f.m,
              wall_height: params_hash['wall_height'].to_f.m,
              sampling_distance: params_hash['sampling_distance'].to_f.m,
              roof_segments: params_hash['roof_segments'].to_i,
              adaptive_sampling: params_hash['adaptive_sampling'] == true,
              smooth_edges: params_hash['smooth_edges'] == true,
              generate_end_caps: params_hash['generate_end_caps'] == true,
              generate_walls: params_hash['generate_walls'] == true,
              generate_roof: params_hash['generate_roof'] == true,
              generate_floor: params_hash['generate_floor'] == true
            }

            puts ">> Parsed Params: #{ruby_params.inspect}"
            callback.call(ruby_params)
          rescue => e
            puts "ERROR: #{e.message}"
            puts e.backtrace.join("\n")
            UI.messagebox("Error: #{e.message}")
          end
        end

        dialog.show
        puts ">> Dialog displayed successfully."
      rescue => e
        puts "ERROR in show_dialog: #{e.message}"
        puts e.backtrace.join("\n")
        UI.messagebox("Error opening dialog: #{e.message}")
      end
    end

    def self.get_html_content
      <<-HTML
      <!DOCTYPE html>
      <html>
      <head>
        <meta charset="utf-8">
        <style>
          body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin: 15px; background-color: #f5f5f5; color: #333; font-size: 13px; }
          .container { background-color: #fff; border: 1px solid #ddd; border-radius: 4px; padding: 15px; }
          h3 { margin-top: 0; color: #444; border-bottom: 1px solid #eee; padding-bottom: 8px; }
          .form-group { margin-bottom: 12px; }
          label { display: block; font-weight: bold; margin-bottom: 4px; }
          input[type="number"] { width: 100%; padding: 6px; box-sizing: border-box; border: 1px solid #ccc; border-radius: 4px; }
          .checkbox-group { margin-top: 8px; }
          .checkbox-label { display: inline-flex; align-items: center; font-weight: normal; margin-right: 15px; margin-bottom: 8px; cursor: pointer; }
          .checkbox-label input { margin-right: 6px; }
          button { width: 100%; padding: 10px; background-color: #0078d4; color: #fff; border: none; border-radius: 4px; font-size: 14px; font-weight: bold; cursor: pointer; margin-top: 10px; }
          button:hover { background-color: #005a9e; }
          .footer { margin-top: 15px; text-align: center; font-size: 11px; color: #666; }
          .status { margin-top: 10px; padding: 8px; background: #f0f0f0; border-radius: 4px; text-align: center; font-size: 12px; color: #666; }
        </style>
      </head>
      <body>
        <div class="container">
          <h3>Tunnel Generation Settings</h3>
          <div class="form-group">
            <label for="roof_radius">Roof Radius (meters)</label>
            <input type="number" id="roof_radius" value="3.0" step="0.1" min="0.1">
          </div>
          <div class="form-group">
            <label for="wall_height">Wall Height (meters)</label>
            <input type="number" id="wall_height" value="2.0" step="0.1" min="0.0">
          </div>
          <div class="form-group">
            <label for="sampling_distance">Sampling Distance (meters)</label>
            <input type="number" id="sampling_distance" value="1.0" step="0.5" min="0.1">
          </div>
          <div class="form-group">
            <label for="roof_segments">Roof Segments</label>
            <input type="number" id="roof_segments" value="12" step="1" min="3">
          </div>
          <div class="checkbox-group">
            <label class="checkbox-label"><input type="checkbox" id="adaptive_sampling" checked> Adaptive Sampling</label><br>
            <label class="checkbox-label"><input type="checkbox" id="smooth_edges" checked> Smooth Edges</label><br>
            <label class="checkbox-label"><input type="checkbox" id="generate_end_caps" checked> Generate End Caps</label><br>
            <label class="checkbox-label"><input type="checkbox" id="generate_walls" checked> Generate Walls</label><br>
            <label class="checkbox-label"><input type="checkbox" id="generate_roof" checked> Generate Roof</label><br>
            <label class="checkbox-label"><input type="checkbox" id="generate_floor" checked> Generate Floor</label>
          </div>
          <button onclick="generateTunnel()">Generate Tunnel</button>
          <div id="status" class="status">Ready</div>
        </div>
        <div class="footer">Tunnel Shell Generator v1.4.0</div>
        <script>
          function generateTunnel() {
            var status = document.getElementById('status');
            status.innerHTML = 'Processing...';
            status.style.color = '#0078d4';

            try {
              var params = {
                roof_radius: parseFloat(document.getElementById('roof_radius').value),
                wall_height: parseFloat(document.getElementById('wall_height').value),
                sampling_distance: parseFloat(document.getElementById('sampling_distance').value),
                roof_segments: parseInt(document.getElementById('roof_segments').value),
                adaptive_sampling: document.getElementById('adaptive_sampling').checked,
                smooth_edges: document.getElementById('smooth_edges').checked,
                generate_end_caps: document.getElementById('generate_end_caps').checked,
                generate_walls: document.getElementById('generate_walls').checked,
                generate_roof: document.getElementById('generate_roof').checked,
                generate_floor: document.getElementById('generate_floor').checked
              };

              window.location = 'skp:generate@' + JSON.stringify(params);
              
              status.innerHTML = '✓ Request sent to SketchUp!';
              status.style.color = '#28a745';
            } catch(e) {
              status.innerHTML = '✗ Error: ' + e.message;
              status.style.color = '#dc3545';
              console.error('Error:', e);
            }
          }
        </script>
      </body>
      </html>
      HTML
    end
  end

  class << self
    def run_generator(params)
      puts ">> run_generator called with params: #{params.inspect}"
      model = Sketchup.active_model
      selection = model.selection

      puts "====================================="
      puts "Tunnel Generator Started"
      puts "====================================="
      puts "1. Selected objects count: #{selection.length}"

      curves = CurveMatching.extract_curves(selection)
      if curves.nil? || curves.length != 2
        UI.messagebox("Error: Please select exactly two curves representing the left and right tunnel edges.")
        puts "ERROR: Curves not found or count is not 2."
        return
      end
      puts "2. Chains found: #{curves.length}"

      c1_edges, c2_edges = curves[0], curves[1]
      puts "3. Edge count - Curve 1: #{c1_edges.length}, Curve 2: #{c2_edges.length}"

      model.start_operation("Generate Tunnel", true)

      begin
        poly1 = CurveMatching.build_polyline(c1_edges)
        poly2 = CurveMatching.build_polyline(c2_edges)
        puts "4. Polylines built. Points: #{poly1.points.length} and #{poly2.points.length}"

        poly1, poly2 = CurveMatching.align_polylines(poly1, poly2)

        samples_t = Sampling.sample_curve(poly1, params)
        puts "5. Sample points count: #{samples_t.length}"

        matched_points = CurveMatching.match_points(poly1, poly2, samples_t)
        left_points = matched_points[:left]
        right_points = matched_points[:right]
        puts "6. Left points: #{left_points.length}, Right points: #{right_points.length}"

        frames = FrameGenerator.generate_frames(left_points, right_points)
        puts "7. Frames generated: #{frames.length}"

        sections = SectionGenerator.generate_sections(frames, params)
        puts "8. Sections generated: #{sections.length}"
        if sections.length > 0
          puts "   - Points per section: #{sections.first.points.length}"
        end

        group = MeshBuilder.build(model, sections, params)
        puts "9. Final group created: #{group.inspect}"

        model.commit_operation
        puts ">> Tunnel generated successfully!"
        UI.messagebox("Tunnel generated successfully!")
      rescue StandardError => e
        model.abort_operation
        UI.messagebox("Error generating tunnel: #{e.message}")
        puts "ERROR: #{e.message}"
        puts e.backtrace.join("\n")
      end
      puts "====================================="
    end

    def register_menu
      unless file_loaded?(__FILE__)
        menu = UI.menu('Extensions')
        menu.add_item('Tunnel Shell Generator') do
          TunnelUI.show_dialog(method(:run_generator))
        end
        file_loaded(__FILE__)
      end
    end
  end

end

TunnelShellGenerator.register_menu
