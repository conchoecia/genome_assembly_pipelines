#!/usr/bin/env python3
"""
Generate interactive HTML report and/or static plots for scaffold classification.

This script creates visualizations showing the relationship between scaffold
mapping percentages to chromosome scaffolds and long-read coverage.
"""

import argparse
import pandas as pd
import sys
import os


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate scaffold classification report and visualizations"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input CSV file with scaffold classification summary"
    )
    parser.add_argument(
        "--alignments",
        help="Optional JSON file with alignment details for dot plots (auto-detected if not provided)"
    )
    parser.add_argument(
        "--output_html",
        help="Output HTML file for interactive report"
    )
    parser.add_argument(
        "--output_static",
        help="Output PNG/PDF file for static plot"
    )
    parser.add_argument(
        "--title",
        default="Scaffold Classification",
        help="Title for the report (default: 'Scaffold Classification')"
    )
    parser.add_argument(
        "--coverage_threshold",
        type=float,
        help="Optional coverage threshold line to display"
    )
    parser.add_argument(
        "--mapping_threshold",
        type=float,
        help="Optional mapping percentage threshold line to display"
    )
    
    return parser.parse_args()


def load_data(input_file):
    """
    Load scaffold classification data from CSV.
    
    Args:
        input_file: Path to input CSV file
        
    Returns:
        pandas DataFrame
    """
    try:
        df = pd.read_csv(input_file)
        print(f"Loaded {len(df)} scaffolds from {input_file}", file=sys.stderr)
        return df
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)


def generate_html_report(df, output_file, title, coverage_threshold=None, mapping_threshold=None, alignment_data=None):
    """
    Generate interactive HTML report using D3.js.
    
    Args:
        df: pandas DataFrame with scaffold data
        output_file: Path to output HTML file
        title: Title for the report
        coverage_threshold: Optional coverage threshold line
        mapping_threshold: Optional mapping percentage threshold line
        alignment_data: Optional dictionary with alignment details for dot plots
    """
    import json
    
    # Prepare data for D3.js
    scaffold_data = []
    for idx, row in df.iterrows():
        # Handle both contig and scaffold naming conventions
        if 'contig_name' in df.columns:
            contig_name = row['contig_name']
            parent_scaffold = row.get('parent_scaffold_name', contig_name)
            display_name = f"{contig_name} ({parent_scaffold})" if contig_name != parent_scaffold else contig_name
            length_col = 'contig_length'
        else:
            contig_name = row['scaffold_name']
            parent_scaffold = None
            display_name = contig_name
            length_col = 'scaffold_length'
        
        scaffold_info = {
            'name': contig_name,
            'display_name': display_name,
            'parent_scaffold': parent_scaffold,
            'length': int(row[length_col]),
            'mapping_pct': float(row['best_match_pct']),
            'coverage': float(row['mean_coverage']),
            'best_chr': str(row['best_match_chr']),
            'total_aligned_pct': float(row.get('percent_aligned', 0)),
            'alignments': alignment_data.get(contig_name, {}) if alignment_data else {}
        }
        scaffold_data.append(scaffold_info)
    
    # Get chromosome columns for heatmap
    chr_cols = [col for col in df.columns if col.endswith('_pct') and col != 'best_match_pct' and col != 'percent_aligned']
    chr_names = [col.replace('_pct', '') for col in chr_cols]
    
    # Prepare heatmap data
    heatmap_data = []
    for idx, row in df.iterrows():
        # Use contig_name or scaffold_name as identifier
        identifier = row.get('contig_name', row.get('scaffold_name'))
        scaffold_heatmap = {
            'scaffold': identifier,
            'values': [float(row[col]) for col in chr_cols]
        }
        heatmap_data.append(scaffold_heatmap)
    
    # Sort heatmap by best match percentage
    identifier_col = 'contig_name' if 'contig_name' in df.columns else 'scaffold_name'
    heatmap_data.sort(key=lambda x: df[df[identifier_col] == x['scaffold']]['best_match_pct'].values[0], reverse=True)
    
    # Calculate summary statistics
    length_col = 'contig_length' if 'contig_length' in df.columns else 'scaffold_length'
    identifier_label = 'Contigs' if 'contig_length' in df.columns else 'Scaffolds'
    total_scaffolds = len(df)
    total_bases = int(df[length_col].sum())
    mean_coverage = float(df['mean_coverage'].mean())
    mean_mapping = float(df['best_match_pct'].mean())
    
    # Generate HTML with D3.js
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Roboto', 'Helvetica', 'Arial', sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
            margin-bottom: 10px;
            font-size: 28px;
        }}
        .subtitle {{
            text-align: center;
            color: #7f8c8d;
            margin-bottom: 30px;
            font-size: 14px;
        }}
        .summary {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 25px;
            border-radius: 8px;
            margin-bottom: 30px;
            color: white;
        }}
        .summary h2 {{
            margin-top: 0;
            font-size: 20px;
            font-weight: 600;
        }}
        .summary-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        .stat-box {{
            text-align: center;
            background: rgba(255, 255, 255, 0.1);
            padding: 15px;
            border-radius: 6px;
        }}
        .stat-value {{
            font-size: 28px;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        .stat-label {{
            font-size: 13px;
            opacity: 0.9;
        }}
        .plot-container {{
            margin: 40px 0;
        }}
        .plot-container h2 {{
            color: #2c3e50;
            font-size: 20px;
            margin-bottom: 20px;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
        }}
        .tooltip {{
            position: absolute;
            padding: 12px;
            background: rgba(0, 0, 0, 0.9);
            color: white;
            border-radius: 6px;
            pointer-events: none;
            font-size: 13px;
            line-height: 1.5;
            display: none;
            max-width: 300px;
            z-index: 1000;
        }}
        .axis-label {{
            font-size: 14px;
            font-weight: 600;
        }}
        svg {{
            display: block;
        }}
        .grid line {{
            stroke: #e0e0e0;
            stroke-opacity: 0.7;
        }}
        .grid path {{
            stroke-width: 0;
        }}
        circle.scatter-point {{
            stroke: white;
            stroke-width: 1.5px;
            cursor: pointer;
            transition: r 0.2s;
        }}
        circle.scatter-point:hover {{
            r: 8 !important;
        }}
        .threshold-line {{
            stroke: #e74c3c;
            stroke-width: 2;
            stroke-dasharray: 5,5;
        }}
        .detail-modal {{
            display: none;
            position: fixed;
            z-index: 2000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            overflow: auto;
            background-color: rgba(0,0,0,0.7);
        }}
        .detail-content {{
            background-color: white;
            margin: 2% auto;
            padding: 20px;
            border-radius: 10px;
            width: 90%;
            max-width: 1200px;
            max-height: 90vh;
            overflow-y: auto;
        }}
        .close-btn {{
            color: #aaa;
            float: right;
            font-size: 28px;
            font-weight: bold;
            cursor: pointer;
            line-height: 20px;
        }}
        .close-btn:hover {{
            color: #000;
        }}
        .detail-section {{
            margin-top: 20px;
        }}
        .detail-section h3 {{
            color: #2c3e50;
            margin-bottom: 10px;
            border-bottom: 2px solid #3498db;
            padding-bottom: 5px;
        }}
        .coverage-bar {{
            fill: #3498db;
            opacity: 0.7;
        }}
        .alignment-block {{
            fill: #e74c3c;
            opacity: 0.6;
            cursor: pointer;
        }}
        .alignment-block:hover {{
            opacity: 0.9;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>
        <div class="subtitle">Competitive Mapping and Coverage Analysis</div>
        
        <div class="summary">
            <h2>Summary Statistics</h2>
            <div class="summary-stats">
                <div class="stat-box">
                    <div class="stat-value">{total_scaffolds:,}</div>
                    <div class="stat-label">Non-Chr {identifier_label}</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{total_bases:,}</div>
                    <div class="stat-label">Total Bases</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{mean_coverage:.2f}x</div>
                    <div class="stat-label">Mean Coverage</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{mean_mapping:.2f}%</div>
                    <div class="stat-label">Mean Mapping %</div>
                </div>
            </div>
        </div>
        
        <div class="plot-container">
            <h2>Mapping vs Coverage Scatter Plot</h2>
            <div id="scatter"></div>
        </div>
        
        <div class="plot-container">
            <h2>Detailed Chromosome Mapping Heatmap</h2>
            <div id="heatmap"></div>
        </div>
    </div>
    
    <div class="tooltip" id="tooltip"></div>
    
    <!-- Detail Modal -->
    <div id="detailModal" class="detail-modal">
        <div class="detail-content">
            <span class="close-btn" onclick="closeDetailModal()">&times;</span>
            <h2 id="modalTitle">Scaffold Details</h2>
            <div id="modalBody"></div>
        </div>
    </div>
    
    <script>
        // Data
        const scaffoldData = {json.dumps(scaffold_data)};
        const heatmapData = {json.dumps(heatmap_data)};
        const chromosomes = {json.dumps(chr_names)};
        const coverageThreshold = {json.dumps(coverage_threshold)};
        const mappingThreshold = {json.dumps(mapping_threshold)};
        
        // Scatter plot
        function createScatterPlot() {{
            const margin = {{top: 20, right: 120, bottom: 60, left: 70}};
            const width = 1000 - margin.left - margin.right;
            const height = 600 - margin.top - margin.bottom;
            
            const svg = d3.select("#scatter")
                .append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", `translate(${{margin.left}},${{margin.top}})`);
            
            // Scales
            const x = d3.scaleLinear()
                .domain([0, 100])
                .range([0, width]);
            
            const maxCoverage = d3.max(scaffoldData, d => d.coverage);
            const y = d3.scaleLinear()
                .domain([0, maxCoverage * 1.1])
                .range([height, 0]);
            
            const colorScale = d3.scaleSequential(d3.interpolateViridis)
                .domain([0, d3.max(scaffoldData, d => d.length)]);
            
            const radiusScale = d3.scaleSqrt()
                .domain([0, d3.max(scaffoldData, d => d.length)])
                .range([3, 10]);
            
            // Grid
            svg.append("g")
                .attr("class", "grid")
                .attr("transform", `translate(0,${{height}})`)
                .call(d3.axisBottom(x).tickSize(-height).tickFormat(""));
            
            svg.append("g")
                .attr("class", "grid")
                .call(d3.axisLeft(y).tickSize(-width).tickFormat(""));
            
            // Axes
            svg.append("g")
                .attr("transform", `translate(0,${{height}})`)
                .call(d3.axisBottom(x))
                .style("font-size", "12px");
            
            svg.append("g")
                .call(d3.axisLeft(y))
                .style("font-size", "12px");
            
            // Axis labels
            svg.append("text")
                .attr("class", "axis-label")
                .attr("text-anchor", "middle")
                .attr("x", width / 2)
                .attr("y", height + 45)
                .text("% of Scaffold Aligned to Best Chr-Scaffold");
            
            svg.append("text")
                .attr("class", "axis-label")
                .attr("text-anchor", "middle")
                .attr("transform", "rotate(-90)")
                .attr("y", -50)
                .attr("x", -height / 2)
                .text("Mean Long-Read Coverage (x)");
            
            // Threshold lines
            if (mappingThreshold !== null) {{
                svg.append("line")
                    .attr("class", "threshold-line")
                    .attr("x1", x(mappingThreshold))
                    .attr("x2", x(mappingThreshold))
                    .attr("y1", 0)
                    .attr("y2", height);
            }}
            
            if (coverageThreshold !== null) {{
                svg.append("line")
                    .attr("class", "threshold-line")
                    .attr("x1", 0)
                    .attr("x2", width)
                    .attr("y1", y(coverageThreshold))
                    .attr("y2", y(coverageThreshold));
            }}
            
            // Tooltip
            const tooltip = d3.select("#tooltip");
            
            // Points
            svg.selectAll("circle")
                .data(scaffoldData)
                .enter()
                .append("circle")
                .attr("class", "scatter-point")
                .attr("cx", d => x(d.mapping_pct))
                .attr("cy", d => y(d.coverage))
                .attr("r", d => radiusScale(d.length))
                .attr("fill", d => colorScale(d.length))
                .attr("opacity", 0.7)
                .on("mouseover", function(event, d) {{
                    d3.select(this).attr("opacity", 1);
                    const hasAlignments = d.alignments && Object.keys(d.alignments).length > 0;
                    tooltip.style("display", "block")
                        .html(`
                            <strong>${{d.display_name || d.name}}</strong><br/>
                            Length: ${{d.length.toLocaleString()}} bp<br/>
                            Best match: ${{d.best_chr}}<br/>
                            Mapping %: ${{d.mapping_pct.toFixed(2)}}%<br/>
                            Coverage: ${{d.coverage.toFixed(2)}}x
                            ${{hasAlignments ? '<br/><em>Click for dot plot</em>' : ''}}
                        `);
                }})
                .on("mousemove", function(event) {{
                    tooltip
                        .style("left", (event.pageX + 15) + "px")
                        .style("top", (event.pageY - 10) + "px");
                }})
                .on("mouseout", function() {{
                    d3.select(this).attr("opacity", 0.7);
                    tooltip.style("display", "none");
                }})
                .on("click", function(event, d) {{
                    if (d.alignments && Object.keys(d.alignments).length > 0) {{
                        showDetailModal(d);
                    }}
                }});
            
            // Color legend
            const legendWidth = 20;
            const legendHeight = 200;
            
            const legendSvg = svg.append("g")
                .attr("transform", `translate(${{width + 20}}, 0)`);
            
            const legendScale = d3.scaleLinear()
                .domain([0, d3.max(scaffoldData, d => d.length)])
                .range([legendHeight, 0]);
            
            const legendAxis = d3.axisRight(legendScale)
                .ticks(5)
                .tickFormat(d => (d / 1000).toFixed(0) + "k");
            
            const defs = svg.append("defs");
            const linearGradient = defs.append("linearGradient")
                .attr("id", "legend-gradient")
                .attr("x1", "0%")
                .attr("x2", "0%")
                .attr("y1", "100%")
                .attr("y2", "0%");
            
            linearGradient.selectAll("stop")
                .data(d3.range(0, 1.01, 0.01))
                .enter()
                .append("stop")
                .attr("offset", d => (d * 100) + "%")
                .attr("stop-color", d => colorScale(d * d3.max(scaffoldData, d => d.length)));
            
            legendSvg.append("rect")
                .attr("width", legendWidth)
                .attr("height", legendHeight)
                .style("fill", "url(#legend-gradient)");
            
            legendSvg.append("g")
                .attr("transform", `translate(${{legendWidth}}, 0)`)
                .call(legendAxis)
                .style("font-size", "11px");
            
            legendSvg.append("text")
                .attr("transform", "rotate(-90)")
                .attr("y", -5)
                .attr("x", -legendHeight / 2)
                .attr("text-anchor", "middle")
                .style("font-size", "11px")
                .text("Length (bp)");
        }}
        
        // Heatmap
        function createHeatmap() {{
            if (heatmapData.length === 0 || chromosomes.length === 0) return;
            
            const margin = {{top: 20, right: 20, bottom: 100, left: 150}};
            const cellHeight = Math.max(15, Math.min(30, 600 / heatmapData.length));
            const cellWidth = Math.max(40, Math.min(80, 800 / chromosomes.length));
            const width = cellWidth * chromosomes.length;
            const height = cellHeight * heatmapData.length;
            
            const svg = d3.select("#heatmap")
                .append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", `translate(${{margin.left}},${{margin.top}})`);
            
            const colorScale = d3.scaleSequential(d3.interpolateYlOrRd)
                .domain([0, 100]);
            
            const tooltip = d3.select("#tooltip");
            
            // X axis (chromosomes)
            const x = d3.scaleBand()
                .domain(chromosomes)
                .range([0, width]);
            
            svg.append("g")
                .attr("transform", `translate(0,${{height}})`)
                .call(d3.axisBottom(x))
                .selectAll("text")
                .attr("transform", "rotate(-45)")
                .style("text-anchor", "end")
                .style("font-size", "11px");
            
            // Y axis (scaffolds) - only show if reasonable number
            const y = d3.scaleBand()
                .domain(heatmapData.map(d => d.scaffold))
                .range([0, height]);
            
            if (heatmapData.length <= 100) {{
                svg.append("g")
                    .call(d3.axisLeft(y))
                    .style("font-size", "10px");
            }}
            
            // Cells
            heatmapData.forEach((scaffold, i) => {{
                scaffold.values.forEach((value, j) => {{
                    svg.append("rect")
                        .attr("x", x(chromosomes[j]))
                        .attr("y", i * cellHeight)
                        .attr("width", cellWidth)
                        .attr("height", cellHeight)
                        .attr("fill", colorScale(value))
                        .attr("stroke", "white")
                        .attr("stroke-width", 0.5)
                        .on("mouseover", function(event) {{
                            d3.select(this).attr("stroke", "black").attr("stroke-width", 2);
                            tooltip.style("display", "block")
                                .html(`
                                    <strong>${{scaffold.scaffold}}</strong><br/>
                                    Chromosome: ${{chromosomes[j]}}<br/>
                                    Mapping: ${{value.toFixed(2)}}%
                                `);
                        }})
                        .on("mousemove", function(event) {{
                            tooltip
                                .style("left", (event.pageX + 15) + "px")
                                .style("top", (event.pageY - 10) + "px");
                        }})
                        .on("mouseout", function() {{
                            d3.select(this).attr("stroke", "white").attr("stroke-width", 0.5);
                            tooltip.style("display", "none");
                        }});
                }});
            }});
        }}
        
        // Modal functions
        function closeDetailModal() {{
            document.getElementById('detailModal').style.display = 'none';
        }}
        
        function showDetailModal(scaffoldData) {{
            const modal = document.getElementById('detailModal');
            const modalTitle = document.getElementById('modalTitle');
            const modalBody = document.getElementById('modalBody');
            
            modalTitle.textContent = `Details for ${{scaffoldData.display_name || scaffoldData.name}}`;
            
            // Build modal content
            let content = `
                <div class="detail-section">
                    <p><strong>Length:</strong> ${{scaffoldData.length.toLocaleString()}} bp</p>
                    <p><strong>Best Match:</strong> ${{scaffoldData.best_chr}} (${{scaffoldData.mapping_pct.toFixed(2)}}% mapped)</p>
                    <p><strong>Coverage:</strong> ${{scaffoldData.coverage.toFixed(2)}}x</p>
                </div>
            `;
            
            // Add dot plot for each chromosome with alignments
            const alignments = scaffoldData.alignments;
            if (alignments && Object.keys(alignments).length > 0) {{
                content += '<div class="detail-section"><h3>Alignment Dot Plots</h3>';
                
                for (const [chrName, chrAlignments] of Object.entries(alignments)) {{
                    if (chrAlignments.length > 0) {{
                        content += `<h4>${{chrName}} (${{chrAlignments.length}} alignment${{chrAlignments.length > 1 ? 's' : ''}})</h4>`;
                        content += `<div id="dotplot_${{chrName.replace(/[^a-zA-Z0-9]/g, '_')}}"></div>`;
                    }}
                }}
                
                content += '</div>';
            }} else {{
                content += '<p><em>No alignment details available</em></p>';
            }}
            
            modalBody.innerHTML = content;
            modal.style.display = 'block';
            
            // Create dot plots after modal is shown
            if (alignments && Object.keys(alignments).length > 0) {{
                for (const [chrName, chrAlignments] of Object.entries(alignments)) {{
                    if (chrAlignments.length > 0) {{
                        const divId = `dotplot_${{chrName.replace(/[^a-zA-Z0-9]/g, '_')}}`;
                        createDotPlot(divId, scaffoldData.name, chrName, chrAlignments, scaffoldData.length);
                    }}
                }}
            }}
        }}
        
        function createDotPlot(divId, scaffoldName, chrName, alignments, scaffoldLength) {{
            const margin = {{top: 20, right: 20, bottom: 50, left: 70}};
            const width = 800 - margin.left - margin.right;
            const height = 300 - margin.top - margin.bottom;
            
            // Clear any existing content
            d3.select(`#${{divId}}`).selectAll("*").remove();
            
            const svg = d3.select(`#${{divId}}`)
                .append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", `translate(${{margin.left}},${{margin.top}})`);
            
            // Find max reference position
            const maxRef = d3.max(alignments, d => d.ref_end);
            
            // Scales
            const xScale = d3.scaleLinear()
                .domain([0, scaffoldLength])
                .range([0, width]);
            
            const yScale = d3.scaleLinear()
                .domain([0, maxRef])
                .range([height, 0]);
            
            // Axes
            const xAxis = d3.axisBottom(xScale)
                .ticks(8)
                .tickFormat(d => (d / 1000).toFixed(0) + "kb");
            
            const yAxis = d3.axisLeft(yScale)
                .ticks(6)
                .tickFormat(d => (d / 1000).toFixed(0) + "kb");
            
            svg.append("g")
                .attr("transform", `translate(0,${{height}})`)
                .call(xAxis);
            
            svg.append("g")
                .call(yAxis);
            
            // Axis labels
            svg.append("text")
                .attr("x", width / 2)
                .attr("y", height + 40)
                .attr("text-anchor", "middle")
                .style("font-size", "12px")
                .text(`${{scaffoldName}} position (bp)`);
            
            svg.append("text")
                .attr("transform", "rotate(-90)")
                .attr("x", -height / 2)
                .attr("y", -50)
                .attr("text-anchor", "middle")
                .style("font-size", "12px")
                .text(`${{chrName}} position (bp)`);
            
            // Tooltip
            const tooltip = d3.select("#tooltip");
            
            // Draw alignment blocks
            alignments.forEach(aln => {{
                // If we have aligned_pairs data, plot those points
                if (aln.aligned_pairs && aln.aligned_pairs.length > 0) {{
                    // Plot aligned pairs as small circles or line segments
                    const color = aln.is_reverse ? "#e74c3c" : "#3498db";
                    
                    // Draw line segments between consecutive aligned pairs
                    for (let i = 0; i < aln.aligned_pairs.length - 1; i++) {{
                        const p1 = aln.aligned_pairs[i];
                        const p2 = aln.aligned_pairs[i + 1];
                        
                        svg.append("line")
                            .attr("x1", xScale(p1[0]))
                            .attr("y1", yScale(p1[1]))
                            .attr("x2", xScale(p2[0]))
                            .attr("y2", yScale(p2[1]))
                            .attr("stroke", color)
                            .attr("stroke-width", 1.5)
                            .attr("opacity", 0.7);
                    }}
                    
                    // Add hover target (invisible wider line)
                    svg.append("line")
                        .attr("x1", xScale(aln.query_start))
                        .attr("x2", xScale(aln.query_end))
                        .attr("y1", yScale(aln.ref_start))
                        .attr("y2", yScale(aln.ref_end))
                        .attr("stroke", "transparent")
                        .attr("stroke-width", 10)
                        .style("cursor", "pointer")
                        .on("mouseover", function(event) {{
                            tooltip.style("display", "block")
                                .html(`
                                    Query: ${{aln.query_start.toLocaleString()}}-${{aln.query_end.toLocaleString()}}<br/>
                                    Ref: ${{aln.ref_start.toLocaleString()}}-${{aln.ref_end.toLocaleString()}}<br/>
                                    MAPQ: ${{aln.mapq}}<br/>
                                    Strand: ${{aln.is_reverse ? 'Reverse' : 'Forward'}}
                                `);
                        }})
                        .on("mousemove", function(event) {{
                            tooltip
                                .style("left", (event.pageX + 15) + "px")
                                .style("top", (event.pageY - 10) + "px");
                        }})
                        .on("mouseout", function() {{
                            tooltip.style("display", "none");
                        }});
                }} else {{
                    // Fallback: draw simple line from start to end
                    svg.append("line")
                        .attr("x1", xScale(aln.query_start))
                        .attr("x2", xScale(aln.query_end))
                        .attr("y1", yScale(aln.ref_start))
                        .attr("y2", yScale(aln.ref_end))
                        .attr("stroke", aln.is_reverse ? "#e74c3c" : "#3498db")
                        .attr("stroke-width", 2)
                        .attr("opacity", 0.6)
                        .on("mouseover", function(event) {{
                            d3.select(this).attr("opacity", 1).attr("stroke-width", 3);
                            tooltip.style("display", "block")
                                .html(`
                                    Query: ${{aln.query_start.toLocaleString()}}-${{aln.query_end.toLocaleString()}}<br/>
                                    Ref: ${{aln.ref_start.toLocaleString()}}-${{aln.ref_end.toLocaleString()}}<br/>
                                    MAPQ: ${{aln.mapq}}<br/>
                                    Strand: ${{aln.is_reverse ? 'Reverse' : 'Forward'}}
                                `);
                        }})
                        .on("mousemove", function(event) {{
                            tooltip
                                .style("left", (event.pageX + 15) + "px")
                                .style("top", (event.pageY - 10) + "px");
                        }})
                        .on("mouseout", function() {{
                            d3.select(this).attr("opacity", 0.6).attr("stroke-width", 2);
                            tooltip.style("display", "none");
                        }});
                }}
            }});
            
            // Legend
            const legend = svg.append("g")
                .attr("transform", `translate(${{width - 100}}, 10)`);
            
            legend.append("line")
                .attr("x1", 0).attr("x2", 20)
                .attr("y1", 0).attr("y2", 0)
                .attr("stroke", "#3498db")
                .attr("stroke-width", 2);
            legend.append("text")
                .attr("x", 25).attr("y", 5)
                .text("Forward")
                .style("font-size", "11px");
            
            legend.append("line")
                .attr("x1", 0).attr("x2", 20)
                .attr("y1", 15).attr("y2", 15)
                .attr("stroke", "#e74c3c")
                .attr("stroke-width", 2);
            legend.append("text")
                .attr("x", 25).attr("y", 20)
                .text("Reverse")
                .style("font-size", "11px");
        }}
        
        // Create visualizations
        createScatterPlot();
        createHeatmap();
    </script>
</body>
</html>"""
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"HTML report written to {output_file}", file=sys.stderr)


def generate_static_plot(df, output_file, title, coverage_threshold=None, mapping_threshold=None):
    """
    Generate static plot using matplotlib.
    
    Args:
        df: pandas DataFrame with scaffold/contig data
        output_file: Path to output PNG/PDF file
        title: Title for the plot
        coverage_threshold: Optional coverage threshold line
        mapping_threshold: Optional mapping percentage threshold line
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("Error: matplotlib is required for static plot generation.", file=sys.stderr)
        print("Install with: pip install matplotlib", file=sys.stderr)
        sys.exit(1)
    
    # Determine column names based on what's available
    length_col = 'contig_length' if 'contig_length' in df.columns else 'scaffold_length'
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Create scatter plot with size and color
    scatter = ax.scatter(
        df['best_match_pct'],
        df['mean_coverage'],
        c=df[length_col],
        s=50,
        alpha=0.6,
        cmap='viridis',
        edgecolors='white',
        linewidth=0.5
    )
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    label_text = 'Contig Length (bp)' if length_col == 'contig_length' else 'Scaffold Length (bp)'
    cbar.set_label(label_text, rotation=270, labelpad=20)
    
    # Add threshold lines if specified
    if mapping_threshold is not None:
        ax.axvline(
            x=mapping_threshold,
            color='red',
            linestyle='--',
            linewidth=2,
            label=f'Mapping threshold: {mapping_threshold}%'
        )
    
    if coverage_threshold is not None:
        ax.axhline(
            y=coverage_threshold,
            color='red',
            linestyle='--',
            linewidth=2,
            label=f'Coverage threshold: {coverage_threshold}x'
        )
    
    # Labels and title
    x_label_text = '% of Contig Aligned to Best Chr-Scaffold' if length_col == 'contig_length' else '% of Scaffold Aligned to Best Chr-Scaffold'
    ax.set_xlabel(x_label_text, fontsize=12)
    ax.set_ylabel('Mean Long-Read Coverage (x)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add legend if thresholds were specified
    if mapping_threshold is not None or coverage_threshold is not None:
        ax.legend()
    
    # Tight layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Static plot written to {output_file}", file=sys.stderr)
    plt.close()


def main():
    """Main function."""
    args = parse_arguments()
    
    # Check that at least one output is specified
    if not args.output_html and not args.output_static:
        print("Error: Must specify at least one of --output_html or --output_static", file=sys.stderr)
        sys.exit(1)
    
    # Load data
    df = load_data(args.input)
    
    # Validate required columns - support both contig and scaffold formats
    if 'contig_name' in df.columns:
        required_cols = ['contig_name', 'contig_length', 'best_match_pct', 'mean_coverage']
    else:
        required_cols = ['scaffold_name', 'scaffold_length', 'best_match_pct', 'mean_coverage']
    
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}", file=sys.stderr)
        sys.exit(1)
    
    # Load alignment details if available
    alignment_data = None
    if args.alignments:
        alignment_file = args.alignments
    else:
        # Auto-detect: look for _alignments.json file next to the input
        # Input is like: path/h12v7_contig_classification_summary.csv
        # Look for: path/h12v7_contig_mapping_stats_alignments.json
        import os
        dir_name = os.path.dirname(args.input)
        base_name = os.path.basename(args.input)
        
        # Extract the assembly name (e.g., h12v7)
        # Assume format: {assembly}_contig_classification_summary.csv
        if '_contig_classification_summary.csv' in base_name:
            assembly_name = base_name.replace('_contig_classification_summary.csv', '')
        elif '_scaffold_classification_summary.csv' in base_name:
            assembly_name = base_name.replace('_scaffold_classification_summary.csv', '')
        else:
            # Fallback: remove .csv
            assembly_name = base_name.replace('.csv', '')
        
        # Construct alignment file path
        alignment_file = os.path.join(dir_name, f'{assembly_name}_contig_mapping_stats_alignments.json')
        print(f"  Looking for alignment file: {alignment_file}", file=sys.stderr)
        
        # Fallback patterns if not found
        if not os.path.exists(alignment_file):
            # Try old scaffold naming for backwards compatibility
            alignment_file = os.path.join(dir_name, f'{assembly_name}_scaffold_mapping_stats_alignments.json')
            print(f"  Trying fallback: {alignment_file}", file=sys.stderr)
        if not os.path.exists(alignment_file):
            alignment_file = os.path.join(dir_name, f'{assembly_name}_mapping_stats_alignments.json')
            print(f"  Trying fallback: {alignment_file}", file=sys.stderr)
        if not os.path.exists(alignment_file):
            alignment_file = os.path.join(dir_name, f'{assembly_name}_alignments.json')
            print(f"  Trying fallback: {alignment_file}", file=sys.stderr)
    
    if os.path.exists(alignment_file):
        print(f"Loading alignment details from {alignment_file}...", file=sys.stderr)
        import json
        try:
            with open(alignment_file, 'r') as f:
                alignment_data = json.load(f)
            print(f"  Loaded alignments for {len(alignment_data)} scaffolds", file=sys.stderr)
        except Exception as e:
            print(f"  Warning: Could not load alignment details: {e}", file=sys.stderr)
    else:
        print(f"  No alignment details found (looked for {alignment_file})", file=sys.stderr)
        print("  Dot plots will not be available in the report", file=sys.stderr)
    
    # Generate HTML report if requested
    if args.output_html:
        generate_html_report(
            df,
            args.output_html,
            args.title,
            args.coverage_threshold,
            args.mapping_threshold,
            alignment_data
        )
    
    # Generate static plot if requested
    if args.output_static:
        generate_static_plot(
            df,
            args.output_static,
            args.title,
            args.coverage_threshold,
            args.mapping_threshold
        )
    
    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
