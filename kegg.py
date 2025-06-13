import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from datetime import datetime
import base64
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import inch
import io
import math
# Remove Streamlit branding and customize appearance
st.set_page_config(
    page_title="AIKEGG Metabolic Pathway Analyzer", 
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Hide Streamlit default elements with CSS
hide_streamlit_style = """
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
.stDeployButton {display:none;}
.stDecoration {display:none;}
[data-testid="stToolbar"] {display: none;}
[data-testid="stHeader"] {display: none;}
[data-testid="stStatusWidget"] {display: none;}
#root > div:nth-child(1) > div > div > div > div > section > div {padding-top: 0rem;}
</style>
"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

@st.cache_data
def get_metabolic_database():
    return {
        'pathways': {
            'hsa04910': {
                'name': 'Insulin signaling pathway',
                'compounds': ['C00031', 'C00267', 'C00089', 'C00103'],
                'compound_names': ['D-Glucose', 'α-D-Glucose', 'Sucrose', 'D-Fructose'],
                'enzymes': ['EC:2.7.1.1', 'EC:2.7.11.1', 'EC:3.1.3.48'],
                'enzyme_names': ['Hexokinase', 'Pyruvate kinase', 'Protein-tyrosine-phosphatase'],
                'reactions': ['R00299', 'R01786', 'R02736'],
                'reaction_names': ['Glucose phosphorylation', 'Pyruvate formation', 'Insulin receptor dephosphorylation'],
                'deltaG': -85.3,
                'disease': 'Diabetes Mellitus',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04910',
                'clinical_significance': 'Critical for glucose homeostasis and insulin sensitivity',
                'affected_organs': ['Pancreas', 'Liver', 'Muscle', 'Adipose tissue'],
                'biomarkers': ['Glucose', 'Insulin', 'C-peptide', 'HbA1c']
            },
            'hsa00010': {
                'name': 'Glycolysis / Gluconeogenesis',
                'compounds': ['C00031', 'C00103', 'C00111', 'C00236'],
                'compound_names': ['D-Glucose', 'D-Fructose', 'Glycerone phosphate', '3-Phospho-D-glycerate'],
                'enzymes': ['EC:2.7.1.1', 'EC:5.3.1.9', 'EC:2.7.1.11'],
                'enzyme_names': ['Hexokinase', 'Glucose-6-phosphate isomerase', '6-phosphofructokinase'],
                'reactions': ['R00299', 'R00771', 'R01068'],
                'reaction_names': ['Glucose phosphorylation', 'Glucose-6-phosphate isomerization', 'Fructose-6-phosphate phosphorylation'],
                'deltaG': -73.3,
                'disease': 'Diabetes Mellitus',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00010',
                'clinical_significance': 'Primary glucose metabolism pathway',
                'affected_organs': ['Liver', 'Muscle', 'Brain', 'Red blood cells'],
                'biomarkers': ['Glucose', 'Lactate', 'Pyruvate']
            },
            'hsa04930': {
                'name': 'Type II diabetes mellitus',
                'compounds': ['C00031', 'C00159', 'C00267', 'C00369'],
                'compound_names': ['D-Glucose', 'D-Mannose', 'α-D-Glucose', 'Starch'],
                'enzymes': ['EC:2.7.1.1', 'EC:3.1.3.48', 'EC:3.1.4.4'],
                'enzyme_names': ['Hexokinase', 'Protein-tyrosine-phosphatase', 'Phospholipase C'],
                'reactions': ['R00299', 'R02736', 'R02740'],
                'reaction_names': ['Glucose phosphorylation', 'Insulin receptor dephosphorylation', 'Phosphoinositide hydrolysis'],
                'deltaG': -68.7,
                'disease': 'Diabetes Mellitus',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04930',
                'clinical_significance': 'Insulin resistance and metabolic dysfunction',
                'affected_organs': ['Pancreas', 'Liver', 'Muscle', 'Adipose tissue'],
                'biomarkers': ['Glucose', 'HbA1c', 'Insulin', 'HOMA-IR']
            },
            'hsa04950': {
                'name': 'Maturity onset diabetes of the young (MODY)',
                'compounds': ['C00031', 'C00103', 'C00267', 'C00668'],
                'compound_names': ['D-Glucose', 'D-Fructose', 'α-D-Glucose', 'ATP'],
                'enzymes': ['EC:2.7.1.1', 'EC:2.7.1.2', 'EC:3.1.3.48'],
                'enzyme_names': ['Hexokinase', 'Glucokinase', 'Protein-tyrosine-phosphatase'],
                'reactions': ['R00299', 'R01786', 'R02736'],
                'reaction_names': ['Glucose phosphorylation', 'Glucose sensing', 'Signal transduction'],
                'deltaG': -78.5,
                'disease': 'Diabetes Mellitus',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04950',
                'clinical_significance': 'Monogenic diabetes with early onset',
                'affected_organs': ['Pancreatic beta cells', 'Liver'],
                'biomarkers': ['Glucose', 'C-peptide', 'Genetic markers']
            },
            'hsa04979': {
                'name': 'Cholesterol metabolism',
                'compounds': ['C00187', 'C00300', 'C00410', 'C01154'],
                'compound_names': ['Cholesterol', 'Creatine', 'Progesterone', 'Cholesteryl ester'],
                'enzymes': ['EC:1.1.1.34', 'EC:1.14.13.70', 'EC:2.3.1.26'],
                'enzyme_names': ['HMG-CoA reductase', 'Cholesterol 7α-hydroxylase', 'Cholesterol acyltransferase'],
                'reactions': ['R00351', 'R01878', 'R02540'],
                'reaction_names': ['Cholesterol biosynthesis', 'Bile acid synthesis', 'Cholesterol esterification'],
                'deltaG': -142.5,
                'disease': 'Dyslipidemia',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04979',
                'clinical_significance': 'Central to lipid homeostasis and cardiovascular risk',
                'affected_organs': ['Liver', 'Intestine', 'Adrenal glands'],
                'biomarkers': ['Total cholesterol', 'LDL-C', 'HDL-C', 'Apo B']
            },
            'hsa03320': {
                'name': 'PPAR signaling pathway',
                'compounds': ['C00187', 'C00410', 'C02530', 'C06427'],
                'compound_names': ['Cholesterol', 'Progesterone', 'Fatty acid', 'Prostaglandin'],
                'enzymes': ['EC:1.13.11.52', 'EC:1.14.19.3', 'EC:1.3.1.22'],
                'enzyme_names': ['Lipoxygenase', 'Fatty acid desaturase', 'Acyl-CoA dehydrogenase'],
                'reactions': ['R01878', 'R02540', 'R04779'],
                'reaction_names': ['Lipid oxidation', 'Fatty acid metabolism', 'Beta-oxidation'],
                'deltaG': -118.9,
                'disease': 'Dyslipidemia',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa03320',
                'clinical_significance': 'Regulates lipid metabolism and inflammation',
                'affected_organs': ['Liver', 'Adipose tissue', 'Muscle'],
                'biomarkers': ['Triglycerides', 'Free fatty acids', 'Inflammatory markers']
            },
            'hsa00561': {
                'name': 'Glycerolipid metabolism',
                'compounds': ['C00044', 'C00116', 'C00162', 'C00422'],
                'compound_names': ['GTP', 'Glycerol', 'Fatty acid', 'Triacylglycerol'],
                'enzymes': ['EC:2.3.1.15', 'EC:3.1.1.3', 'EC:3.1.1.4'],
                'enzyme_names': ['Glycerol-3-phosphate acyltransferase', 'Triacylglycerol lipase', 'Phospholipase A2'],
                'reactions': ['R00256', 'R01046', 'R01049'],
                'reaction_names': ['Glycerol phosphorylation', 'Triacylglycerol synthesis', 'Lipolysis'],
                'deltaG': -92.4,
                'disease': 'Dyslipidemia',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00561',
                'clinical_significance': 'Triglyceride synthesis and breakdown',
                'affected_organs': ['Liver', 'Adipose tissue', 'Muscle'],
                'biomarkers': ['Triglycerides', 'Glycerol', 'Free fatty acids']
            },
            'hsa04960': {
                'name': 'Aldosterone-regulated sodium reabsorption',
                'compounds': ['C00009', 'C00080', 'C01330', 'C00788'],
                'compound_names': ['Orthophosphate', 'H+', 'Sodium', 'L-Histidine'],
                'enzymes': ['EC:1.6.3.1', 'EC:3.6.3.9', 'EC:4.2.1.1'],
                'enzyme_names': ['NAD(P)H oxidase', 'Na+/K+-ATPase', 'Carbonic anhydrase'],
                'reactions': ['R00114', 'R00115', 'R00116'],
                'reaction_names': ['Sodium transport', 'Potassium exchange', 'pH regulation'],
                'deltaG': -45.2,
                'disease': 'Chronic Kidney Disease',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04960',
                'clinical_significance': 'Electrolyte balance and blood pressure regulation',
                'affected_organs': ['Kidney', 'Adrenal cortex'],
                'biomarkers': ['Sodium', 'Potassium', 'Aldosterone', 'Renin']
            },
            'hsa04961': {
                'name': 'Endocrine and other factor-regulated calcium reabsorption',
                'compounds': ['C00076', 'C01353', 'C00005', 'C00009'],
                'compound_names': ['Calcium', 'Calcitriol', 'NADPH', 'Orthophosphate'],
                'enzymes': ['EC:3.6.3.8', 'EC:1.1.1.37', 'EC:4.2.1.1'],
                'enzyme_names': ['Ca2+-ATPase', 'Malate dehydrogenase', 'Carbonic anhydrase'],
                'reactions': ['R00114', 'R00115', 'R00116'],
                'reaction_names': ['Calcium transport', 'Energy metabolism', 'Acid-base balance'],
                'deltaG': -38.7,
                'disease': 'Chronic Kidney Disease',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04961',
                'clinical_significance': 'Calcium homeostasis and bone metabolism',
                'affected_organs': ['Kidney', 'Parathyroid', 'Bone'],
                'biomarkers': ['Calcium', 'Phosphate', 'PTH', 'Vitamin D']
            },
            'hsa04966': {
                'name': 'Collecting duct acid secretion',
                'compounds': ['C00014', 'C00058', 'C00238', 'C01353'],
                'compound_names': ['Ammonia', 'Formate', 'Potassium', 'Calcitriol'],
                'enzymes': ['EC:3.6.3.1', 'EC:4.2.1.1', 'EC:1.1.1.37'],
                'enzyme_names': ['H+-ATPase', 'Carbonic anhydrase', 'Malate dehydrogenase'],
                'reactions': ['R00114', 'R00115', 'R00116'],
                'reaction_names': ['Proton secretion', 'Bicarbonate formation', 'Acid-base regulation'],
                'deltaG': -52.1,
                'disease': 'Chronic Kidney Disease',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04966',
                'clinical_significance': 'Renal acid-base balance and pH homeostasis',
                'affected_organs': ['Kidney collecting duct'],
                'biomarkers': ['pH', 'Bicarbonate', 'Anion gap', 'Ammonia']
            }
        }
    }

@st.cache_data
def get_reference_ranges():
    return {
        'glucose': {'min': 70, 'max': 100, 'unit': 'mg/dL'},
        'bun': {'min': 6, 'max': 24, 'unit': 'mg/dL'},
        'creatinine': {'min': 0.7, 'max': 1.3, 'unit': 'mg/dL'},
        'sodium': {'min': 136, 'max': 145, 'unit': 'mEq/L'},
        'potassium': {'min': 3.5, 'max': 5.1, 'unit': 'mEq/L'},
        'cholesterol': {'min': 0, 'max': 200, 'unit': 'mg/dL'},
        'hdl': {'min': 40, 'max': 200, 'unit': 'mg/dL'},
        'ldl': {'min': 0, 'max': 100, 'unit': 'mg/dL'},
        'vldl': {'min': 5, 'max': 40, 'unit': 'mg/dL'},
        'triglycerides': {'min': 0, 'max': 150, 'unit': 'mg/dL'},
        'alt': {'min': 7, 'max': 56, 'unit': 'U/L'},
        'ast': {'min': 10, 'max': 40, 'unit': 'U/L'},
        'alp': {'min': 44, 'max': 147, 'unit': 'IU/L'},
        'bilirubin': {'min': 0.1, 'max': 1.2, 'unit': 'mg/dL'},
        'albumin': {'min': 3.5, 'max': 5.0, 'unit': 'g/dL'},
        'uric_acid': {'min': 3.4, 'max': 7.0, 'unit': 'mg/dL'},
        'calcium': {'min': 8.5, 'max': 10.5, 'unit': 'mg/dL'},
        'phosphate': {'min': 2.5, 'max': 4.5, 'unit': 'mg/dL'},
        'hba1c': {'min': 4.0, 'max': 5.6, 'unit': '%'},
        'insulin': {'min': 2.6, 'max': 24.9, 'unit': 'μIU/mL'},
        'c_peptide': {'min': 1.1, 'max': 4.4, 'unit': 'ng/mL'}
    }

def calculate_z_score(value, mean, std):
    return (value - mean) / std

def get_pathway_perturbation_score(lab_values, metabolic_db, ref_ranges):
    scores = {}
    
    for pathway_id, pathway in metabolic_db['pathways'].items():
        perturbation_score = 0
        relevant_values = 0
        
        # Glucose-related pathways
        if 'glucose' in lab_values and 'C00031' in pathway['compounds']:
            reference = ref_ranges['glucose']
            mean = (reference['min'] + reference['max']) / 2
            std = (reference['max'] - reference['min']) / 4
            perturbation_score += abs(calculate_z_score(lab_values['glucose'], mean, std))
            relevant_values += 1
        
        # Cholesterol-related pathways
        if 'cholesterol' in lab_values and 'C00187' in pathway['compounds']:
            reference = ref_ranges['cholesterol']
            mean = (reference['min'] + reference['max']) / 2
            std = (reference['max'] - reference['min']) / 4
            perturbation_score += abs(calculate_z_score(lab_values['cholesterol'], mean, std))
            relevant_values += 1
        
        # Triglycerides-related pathways
        if 'triglycerides' in lab_values and 'C00162' in pathway['compounds']:
            reference = ref_ranges['triglycerides']
            mean = (reference['min'] + reference['max']) / 2
            std = (reference['max'] - reference['min']) / 4
            perturbation_score += abs(calculate_z_score(lab_values['triglycerides'], mean, std))
            relevant_values += 1
        
        # Kidney function markers
        if 'creatinine' in lab_values and pathway['disease'] == 'Chronic Kidney Disease':
            reference = ref_ranges['creatinine']
            mean = (reference['min'] + reference['max']) / 2
            std = (reference['max'] - reference['min']) / 4
            perturbation_score += abs(calculate_z_score(lab_values['creatinine'], mean, std))
            relevant_values += 1
            
        if 'bun' in lab_values and pathway['disease'] == 'Chronic Kidney Disease':
            reference = ref_ranges['bun']
            mean = (reference['min'] + reference['max']) / 2
            std = (reference['max'] - reference['min']) / 4
            perturbation_score += abs(calculate_z_score(lab_values['bun'], mean, std))
            relevant_values += 1
        
        # Additional biomarkers
        if 'calcium' in lab_values and 'C00076' in pathway['compounds']:
            reference = ref_ranges['calcium']
            mean = (reference['min'] + reference['max']) / 2
            std = (reference['max'] - reference['min']) / 4
            perturbation_score += abs(calculate_z_score(lab_values['calcium'], mean, std))
            relevant_values += 1
            
        if 'phosphate' in lab_values and 'C00009' in pathway['compounds']:
            reference = ref_ranges['phosphate']
            mean = (reference['min'] + reference['max']) / 2
            std = (reference['max'] - reference['min']) / 4
            perturbation_score += abs(calculate_z_score(lab_values['phosphate'], mean, std))
            relevant_values += 1
        
        # Liver function markers
        if 'alt' in lab_values and pathway['disease'] in ['Diabetes Mellitus', 'Dyslipidemia']:
            reference = ref_ranges['alt']
            mean = (reference['min'] + reference['max']) / 2
            std = (reference['max'] - reference['min']) / 4
            perturbation_score += abs(calculate_z_score(lab_values['alt'], mean, std))
            relevant_values += 1
            
        if 'ast' in lab_values and pathway['disease'] in ['Diabetes Mellitus', 'Dyslipidemia']:
            reference = ref_ranges['ast']
            mean = (reference['min'] + reference['max']) / 2
            std = (reference['max'] - reference['min']) / 4
            perturbation_score += abs(calculate_z_score(lab_values['ast'], mean, std))
            relevant_values += 1
        
        # Calculate thermodynamic score
        thermodynamic_score = abs(pathway['deltaG']) / 100
        
        # Store scores in dictionary
        scores[pathway_id] = {
            'perturbation_score': perturbation_score / max(relevant_values, 1),
            'thermodynamic_score': thermodynamic_score,
            'total_score': (perturbation_score / max(relevant_values, 1)) + thermodynamic_score,
            'pathway': pathway,
            'relevant_markers': relevant_values
        }
    
    return scores

def get_affected_pathways_analysis(lab_values, metabolic_db, ref_ranges):
    """Analyze which pathways are directly affected and at risk"""
    directly_affected = []
    at_risk = []
    
    pathway_scores = get_pathway_perturbation_score(lab_values, metabolic_db, ref_ranges)
    
    for pathway_id, data in pathway_scores.items():
        if data['perturbation_score'] > 1.5:  # High perturbation
            directly_affected.append({
                'pathway_id': pathway_id,
                'name': data['pathway']['name'],
                'disease': data['pathway']['disease'],
                'score': data['total_score'],
                'affected_biomarkers': get_affected_biomarkers(lab_values, data['pathway'], ref_ranges)
            })
        elif data['perturbation_score'] > 0.5:  # Moderate perturbation
            at_risk.append({
                'pathway_id': pathway_id,
                'name': data['pathway']['name'],
                'disease': data['pathway']['disease'],
                'score': data['total_score'],
                'risk_factors': get_risk_factors(lab_values, data['pathway'], ref_ranges)
            })
    
    return directly_affected, at_risk

def get_affected_biomarkers(lab_values, pathway, ref_ranges):
    """Get biomarkers that are abnormal and related to the pathway"""
    affected = []
    for biomarker in pathway['biomarkers']:
        biomarker_key = biomarker.lower().replace('-', '_').replace(' ', '_')
        if biomarker_key in lab_values and biomarker_key in ref_ranges:
            value = lab_values[biomarker_key]
            ref = ref_ranges[biomarker_key]
            if value < ref['min'] or value > ref['max']:
                status = 'Low' if value < ref['min'] else 'High'
                affected.append(f"{biomarker}: {status}")
    return affected

def get_risk_factors(lab_values, pathway, ref_ranges):
    """Get risk factors for pathways at risk"""
    risk_factors = []
    for biomarker in pathway['biomarkers']:
        biomarker_key = biomarker.lower().replace('-', '_').replace(' ', '_')
        if biomarker_key in lab_values and biomarker_key in ref_ranges:
            value = lab_values[biomarker_key]
            ref = ref_ranges[biomarker_key]
            # Check if value is in upper 75% of normal range (potential risk)
            upper_75 = ref['min'] + 0.75 * (ref['max'] - ref['min'])
            if ref['min'] <= value <= ref['max'] and value > upper_75:
                risk_factors.append(f"{biomarker}: Upper normal range")
    return risk_factors

def create_overview_chart(sorted_pathways):
    pathway_names = [data[1]['pathway']['name'][:20] + '...' if len(data[1]['pathway']['name']) > 20 else data[1]['pathway']['name'] for data in sorted_pathways[:5]]
    pathway_scores = [data[1]['total_score'] for data in sorted_pathways[:5]]
    
    colors = ['#EF4444', '#F56565', '#FB923C', '#22C55E', '#3B82F6']
    
    fig = go.Figure(data=[
        go.Bar(
            y=pathway_names,
            x=pathway_scores,
            orientation='h',
            marker_color=colors,
            text=[f"{score:.2f}" for score in pathway_scores],
            textposition='outside'
        )
    ])
    
    fig.update_layout(
        title="Top 5 Affected Metabolic Pathways",
        xaxis_title="Perturbation Score",
        height=400,
        margin=dict(l=200)
    )
    
    return fig

def create_thermodynamics_chart(sorted_pathways):
    pathway_names = [data[1]['pathway']['name'][:15] + '...' if len(data[1]['pathway']['name']) > 15 else data[1]['pathway']['name'] for data in sorted_pathways]
    delta_g = [data[1]['pathway']['deltaG'] for data in sorted_pathways]
    
    colors = ['#48BB78' if dg < -50 else '#F56565' if dg > -20 else '#718096' for dg in delta_g]
    
    fig = go.Figure(data=[
        go.Bar(
            x=pathway_names,
            y=delta_g,
            marker_color=colors,
            text=[f"{dg:.1f}" for dg in delta_g],
            textposition='outside'
        )
    ])
    
    fig.update_layout(
        title="Thermodynamic Favorability of Metabolic Pathways",
        yaxis_title="Gibbs Free Energy (kJ/mol)",
        height=400
    )
    
    return fig

def create_flux_chart(sorted_pathways):
    pathway_names = [data[1]['pathway']['name'][:12] + '...' if len(data[1]['pathway']['name']) > 12 else data[1]['pathway']['name'] for data in sorted_pathways[:5]]
    
    normal_flux = [100] * 5
    perturbed_flux = [100 * (1 + (data[1]['perturbation_score'] * 0.3)) for data in sorted_pathways[:5]]
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatterpolar(
        r=normal_flux,
        theta=pathway_names,
        fill='toself',
        name='Normal Flux',
        line_color='rgba(72, 187, 120, 1)'
    ))
    
    fig.add_trace(go.Scatterpolar(
        r=perturbed_flux,
        theta=pathway_names,
        fill='toself',
        name='Perturbed Flux',
        line_color='rgba(245, 101, 101, 1)'
    ))
    
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 150]
            )),
        showlegend=True,
        title="Metabolic Flux Analysis - Normal vs Perturbed States",
        height=500
    )
    
    return fig

def create_kegg_style_pathway_chart(pathway_data, selected_pathway_id=None):
    """Create a KEGG-style metabolic pathway visualization"""
    
    if selected_pathway_id and selected_pathway_id in pathway_data:
        # Show detailed view of single pathway
        pathway = pathway_data[selected_pathway_id]['pathway']
        
        # Create a more detailed single pathway view
        fig = go.Figure()
        
        # Define positions for compounds, enzymes, and reactions
        compounds = pathway['compounds'][:4]
        compound_names = pathway['compound_names'][:4]
        enzymes = pathway['enzymes'][:3]
        enzyme_names = pathway['enzyme_names'][:3]
        reactions = pathway['reactions'][:3]
        reaction_names = pathway['reaction_names'][:3]
        
        # Positions for KEGG-like layout
        compound_positions = [(0, 2), (2, 3), (4, 2), (6, 1)]
        enzyme_positions = [(1, 1), (3, 1), (5, 1)]
        reaction_positions = [(1, 2.5), (3, 2.5), (5, 2.5)]
        
        # Add compounds (circles)
        for i, (comp_id, comp_name, pos) in enumerate(zip(compounds, compound_names, compound_positions)):
            fig.add_trace(go.Scatter(
                x=[pos[0]], y=[pos[1]],
                mode='markers+text',
                marker=dict(size=40, color='lightblue', line=dict(width=2, color='darkblue')),
                text=comp_name.split()[0],  # Short name
                textposition="middle center",
                name=f"Compound",
                hovertext=f"{comp_name} ({comp_id})",
                showlegend=False
            ))
        
        # Add enzymes (rectangles - represented as squares)
        for i, (enz_id, enz_name, pos) in enumerate(zip(enzymes, enzyme_names, enzyme_positions)):
            fig.add_trace(go.Scatter(
                x=[pos[0]], y=[pos[1]],
                mode='markers+text',
                marker=dict(size=35, color='lightgreen', symbol='square', line=dict(width=2, color='darkgreen')),
                text=f"E{i+1}",
                textposition="middle center",
                name=f"Enzyme",
                hovertext=f"{enz_name} ({enz_id})",
                showlegend=False
            ))
        
        # Add reaction arrows
        for i, (pos1, pos2) in enumerate(zip(compound_positions[:-1], compound_positions[1:])):
            # Arrow from compound to compound via enzyme
            enz_pos = enzyme_positions[i] if i < len(enzyme_positions) else enzyme_positions[-1]
            
            # Line from compound to enzyme
            fig.add_trace(go.Scatter(
                x=[pos1[0], enz_pos[0]], y=[pos1[1], enz_pos[1]],
                mode='lines',
                line=dict(width=2, color='red'),
                showlegend=False,
                hoverinfo='skip'
            ))
            
            # Line from enzyme to next compound
            if i+1 < len(compound_positions):
                fig.add_trace(go.Scatter(
                    x=[enz_pos[0], pos2[0]], y=[enz_pos[1], pos2[1]],
                    mode='lines',
                    line=dict(width=2, color='red'),
                    showlegend=False,
                    hoverinfo='skip'
                ))
        
        # Add pathway title and info
        fig.update_layout(
            title=f"KEGG-Style View: {pathway['name']}",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.5, 6.5]),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[0, 4]),
            height=500,
            annotations=[
                dict(text=f"Disease: {pathway['disease']}", x=0, y=3.5, showarrow=False, font=dict(size=12)),
                dict(text=f"ΔG: {pathway['deltaG']} kJ/mol", x=0, y=3.2, showarrow=False, font=dict(size=12)),
                dict(text="🔵 Compounds  🟩 Enzymes  ➡️ Reactions", x=3, y=0.2, showarrow=False, font=dict(size=10))
            ]
        )
        
    else:
        # Overview network showing pathway connections
        fig = create_pathway_overview_network(pathway_data)
    
    return fig

def create_pathway_overview_network(pathway_data):
    """Create an overview network of pathway relationships"""
    import networkx as nx
    
    # Group pathways by disease
    disease_groups = {}
    for pathway_id, data in pathway_data.items():
        disease = data['pathway']['disease']
        if disease not in disease_groups:
            disease_groups[disease] = []
        disease_groups[disease].append((pathway_id, data))
    
    fig = go.Figure()
    
    # Color mapping for diseases
    colors = {
        'Diabetes Mellitus': '#FF6B6B',
        'Dyslipidemia': '#4ECDC4', 
        'Chronic Kidney Disease': '#45B7D1'
    }
    
    # Position diseases in clusters
    disease_positions = {
        'Diabetes Mellitus': (0, 0),
        'Dyslipidemia': (3, 0),
        'Chronic Kidney Disease': (1.5, 2.5)
    }
    
    # Add disease clusters
    for disease, pathways in disease_groups.items():
        center_pos = disease_positions.get(disease, (0, 0))
        
        # Add disease label
        fig.add_trace(go.Scatter(
            x=[center_pos[0]], y=[center_pos[1] + 0.8],
            mode='text',
            text=disease,
            textfont=dict(size=16, color=colors.get(disease, '#333')),
            showlegend=False,
            hoverinfo='skip'
        ))
        
        # Add pathways in cluster
        for i, (pathway_id, data) in enumerate(pathways):
            angle = i * (2 * 3.14159 / len(pathways))
            radius = 0.5
            x = center_pos[0] + radius * math.cos(angle)
            y = center_pos[1] + radius * math.sin(angle)
            
            fig.add_trace(go.Scatter(
                x=[x], y=[y],
                mode='markers+text',
                marker=dict(
                    size=max(20, min(40, data['total_score'] * 10)),
                    color=colors.get(disease, '#333'),
                    opacity=0.7,
                    line=dict(width=2, color='white')
                ),
                text=data['pathway']['name'][:8] + '...' if len(data['pathway']['name']) > 8 else data['pathway']['name'],
                textposition="bottom center",
                textfont=dict(size=8),
                name=disease,
                hovertext=f"{data['pathway']['name']}<br>Score: {data['total_score']:.2f}<br>Click for detailed view",
                showlegend=i == 0  # Only show legend for first pathway of each disease
            ))
    
    # Add connections between related pathways
    for disease1, pathways1 in disease_groups.items():
        for disease2, pathways2 in disease_groups.items():
            if disease1 != disease2:
                # Add connection line between disease centers
                pos1 = disease_positions[disease1]
                pos2 = disease_positions[disease2]
                fig.add_trace(go.Scatter(
                    x=[pos1[0], pos2[0]], y=[pos1[1], pos2[1]],
                    mode='lines',
                    line=dict(width=1, color='lightgray', dash='dash'),
                    showlegend=False,
                    hoverinfo='skip'
                ))
    
    fig.update_layout(
        title="KEGG-Style Metabolic Pathway Network Overview",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-1, 4]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-1, 4]),
        height=600,
        hovermode='closest'
    )
    
    return fig

def display_pathway_details(pathway_data, pathway_id):
    """Display detailed information about a specific pathway"""
    data = pathway_data[pathway_id]
    pathway = data['pathway']
    
    st.subheader(f"🔬 {pathway['name']}")
    st.write(f"**Associated Disease:** {pathway['disease']}")
    st.write(f"**Clinical Significance:** {pathway['clinical_significance']}")
    st.write(f"**KEGG Map:** [View Pathway]({pathway['kegg_map_url']})")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**Key Compounds (4):**")
        for i, (compound_id, compound_name) in enumerate(zip(pathway['compounds'][:4], pathway['compound_names'][:4])):
            st.write(f"• {compound_name} ({compound_id})")
        
        st.write("**Enzymes (3):**")
        for i, (enzyme_id, enzyme_name) in enumerate(zip(pathway['enzymes'][:3], pathway['enzyme_names'][:3])):
            st.write(f"• {enzyme_name} ({enzyme_id})")
    
    with col2:
        st.write("**Reactions (3):**")
        for i, (reaction_id, reaction_name) in enumerate(zip(pathway['reactions'][:3], pathway['reaction_names'][:3])):
            st.write(f"• {reaction_name} ({reaction_id})")
        
        st.write("**Affected Organs:**")
        for organ in pathway['affected_organs']:
            st.write(f"• {organ}")
    
    st.write("**Relevant Biomarkers:**")
    biomarker_text = ", ".join(pathway['biomarkers'])
    st.write(biomarker_text)
    
    # Thermodynamic information
    favorability = "Highly Favorable" if pathway['deltaG'] < -50 else "Moderately Favorable" if pathway['deltaG'] < -20 else "Slightly Favorable" if pathway['deltaG'] < 0 else "Unfavorable"
    st.write(f"**Thermodynamic Status:** ΔG = {pathway['deltaG']:.1f} kJ/mol ({favorability})")

def generate_enhanced_report_content(lab_values, sorted_pathways, patient_info, ref_ranges):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Analyze directly affected and at-risk pathways
    directly_affected, at_risk = get_affected_pathways_analysis(lab_values, get_metabolic_database(), ref_ranges)
    
    lab_values_text = ""
    for param, value in lab_values.items():
        ref = ref_ranges.get(param, {})
        status = "Normal"
        if ref:
            if value < ref['min']:
                status = "Low ⬇️"
            elif value > ref['max']:
                status = "High ⬆️"
        
        lab_values_text += f"{param.upper()}: {value} {ref.get('unit', '')} (Ref: {ref.get('min', 'N/A')}-{ref.get('max', 'N/A')}) - {status}\n"
    
    # Directly affected pathways analysis
    directly_affected_text = ""
    if directly_affected:
        directly_affected_text = "DIRECTLY AFFECTED PATHWAYS:\n"
        for pathway in directly_affected:
            directly_affected_text += f"• {pathway['name']} ({pathway['disease']}) - Score: {pathway['score']:.2f}\n"
            if pathway['affected_biomarkers']:
                directly_affected_text += f"  Abnormal biomarkers: {', '.join(pathway['affected_biomarkers'])}\n"
    else:
        directly_affected_text = "DIRECTLY AFFECTED PATHWAYS: None detected\n"
    
    # At-risk pathways analysis
    at_risk_text = ""
    if at_risk:
        at_risk_text = "\nPATHWAYS AT RISK:\n"
        for pathway in at_risk:
            at_risk_text += f"• {pathway['name']} ({pathway['disease']}) - Score: {pathway['score']:.2f}\n"
            if pathway['risk_factors']:
                at_risk_text += f"  Risk factors: {', '.join(pathway['risk_factors'])}\n"
    else:
        at_risk_text = "\nPATHWAYS AT RISK: None detected\n"
    
    # Top pathway details
    top_pathway = sorted_pathways[0][1]['pathway']
    pathway_details = f"""
TOP PRIORITY PATHWAY DETAILS:
Name: {top_pathway['name']}
Disease Association: {top_pathway['disease']}
Key Compounds: {', '.join(top_pathway['compound_names'][:4])}
Key Enzymes: {', '.join(top_pathway['enzyme_names'][:3])}
Key Reactions: {', '.join(top_pathway['reaction_names'][:3])}
Thermodynamic Status: ΔG = {top_pathway['deltaG']:.1f} kJ/mol
KEGG Map: {top_pathway['kegg_map_url']}
"""
    
    abnormal_count = sum(1 for param, value in lab_values.items() 
                        if param in ref_ranges and 
                        (value < ref_ranges[param]['min'] or value > ref_ranges[param]['max']))
    
    report_content = f"""AIKEGG METABOLIC PATHWAY ANALYSIS REPORT
Generated: {timestamp}

PATIENT INFORMATION
Patient ID: {patient_info.get('patient_id', 'Not specified')}
Age: {patient_info.get('age', 'Not specified')}
Gender: {patient_info.get('gender', 'Not specified')}

LABORATORY VALUES
{lab_values_text}

PATHWAY PERTURBATION ANALYSIS
{directly_affected_text}
{at_risk_text}

{pathway_details}

CLINICAL RECOMMENDATIONS
• Immediate attention required for: {len(directly_affected)} pathway(s)
• Monitoring recommended for: {len(at_risk)} pathway(s)
• Consider follow-up testing for biomarkers in affected pathways

SUMMARY
Total Parameters Analyzed: {len(lab_values)}
Abnormal Values Detected: {abnormal_count}
Pathways Directly Affected: {len(directly_affected)}
Pathways At Risk: {len(at_risk)}
"""
    
    return report_content
def generate_pdf_report(lab_values, sorted_pathways, patient_info, ref_ranges):
    """Generate PDF report using ReportLab"""
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    styles = getSampleStyleSheet()
    story = []
    
    # Title
    title_style = ParagraphStyle('CustomTitle', parent=styles['Heading1'], fontSize=18, spaceAfter=30)
    story.append(Paragraph("AIKEGG METABOLIC PATHWAY ANALYSIS REPORT", title_style))
    story.append(Spacer(1, 12))
    
    # Patient Info
    story.append(Paragraph("PATIENT INFORMATION", styles['Heading2']))
    patient_data = [
        ['Patient ID:', patient_info.get('patient_id', 'Not specified')],
        ['Age:', str(patient_info.get('age', 'Not specified'))],
        ['Gender:', patient_info.get('gender', 'Not specified')],
        ['Report Date:', datetime.now().strftime("%Y-%m-%d %H:%M:%S")]
    ]
    patient_table = Table(patient_data, colWidths=[2*inch, 3*inch])
    patient_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), colors.lightgrey),
        ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 12),
    ]))
    story.append(patient_table)
    story.append(Spacer(1, 20))
    
    # Lab Values
    story.append(Paragraph("LABORATORY VALUES", styles['Heading2']))
    lab_data = [['Parameter', 'Value', 'Reference Range', 'Status']]
    
    for param, value in lab_values.items():
        ref = ref_ranges.get(param, {})
        status = "Normal"
        if ref:
            if value < ref['min']:
                status = "Low"
            elif value > ref['max']:
                status = "High"
        
        lab_data.append([
            param.upper(),
            f"{value} {ref.get('unit', '')}",
            f"{ref.get('min', 'N/A')}-{ref.get('max', 'N/A')} {ref.get('unit', '')}",
            status
        ])
    
    lab_table = Table(lab_data, colWidths=[1.5*inch, 1.5*inch, 1.5*inch, 1*inch])
    lab_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black)
    ]))
    story.append(lab_table)
    story.append(Spacer(1, 20))
    
    # Top Pathways
    story.append(Paragraph("TOP AFFECTED PATHWAYS", styles['Heading2']))
    pathway_data = [['Rank', 'Pathway Name', 'Disease', 'Score']]
    
    for i, (pathway_id, data) in enumerate(sorted_pathways[:5], 1):
        pathway_data.append([
            str(i),
            data['pathway']['name'],
            data['pathway']['disease'],
            f"{data['total_score']:.2f}"
        ])
    
    pathway_table = Table(pathway_data, colWidths=[0.5*inch, 3*inch, 1.5*inch, 1*inch])
    pathway_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black)
    ]))
    story.append(pathway_table)
    
    doc.build(story)
    buffer.seek(0)
    return buffer
def create_pathway_table(sorted_pathways, directly_affected, at_risk):
    """Create an Excel-like table for pathway analysis"""
    
    # Prepare data for the table
    table_data = []
    
    for i, (pathway_id, data) in enumerate(sorted_pathways, 1):
        # Determine status
        status = "Normal"
        status_color = "🟢"
        
        # Check if pathway is directly affected
        for affected in directly_affected:
            if affected['pathway_id'] == pathway_id:
                status = "Directly Affected"
                status_color = "🔴"
                break
        
        # Check if pathway is at risk
        if status == "Normal":
            for risk in at_risk:
                if risk['pathway_id'] == pathway_id:
                    status = "At Risk"
                    status_color = "🟡"
                    break
        
        table_data.append({
            'Rank': i,
            'Status': f"{status_color} {status}",
            'Pathway Name': data['pathway']['name'],
            'Disease Association': data['pathway']['disease'],
            'Perturbation Score': f"{data['perturbation_score']:.2f}",
            'Thermodynamic Score': f"{data['thermodynamic_score']:.2f}",
            'Total Score': f"{data['total_score']:.2f}",
            'Affected Organs': ', '.join(data['pathway']['affected_organs'][:3]),
            'Key Biomarkers': ', '.join(data['pathway']['biomarkers'][:3]),
            'ΔG (kJ/mol)': f"{data['pathway']['deltaG']:.1f}"
        })
    
    return table_data
# Update the main function with enhanced pathway analysis
def main():
    st.title("🧬 AIKEGG - AI-Enhanced Metabolic Pathway Analysis")
    st.markdown("Advanced AI platform that maps patient biochemical data to KEGG metabolic pathways for precision diagnostics")
    
    metabolic_db = get_metabolic_database()
    ref_ranges = get_reference_ranges()
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.header("Patient Information")
        patient_id = st.text_input("Patient ID", placeholder="Enter patient identifier")
        age = st.number_input("Age", min_value=1, max_value=120, value=None)
        gender = st.selectbox("Gender", ["", "Male", "Female", "Other"])
        
        st.header("Disease Modules")
        st.info("Focus on Diabetes, Dyslipidemia, and CKD pathways based on clinical prevalence")
        
        # Metabolic Panel
        with st.expander("🩸 Metabolic Panel", expanded=True):
            glucose = st.number_input("Fasting Glucose (mg/dL)", min_value=0.0, value=None, step=0.1)
            creatinine = st.number_input("Serum Creatinine (mg/dL)", min_value=0.0, value=None, step=0.01)
            bun = st.number_input("Blood Urea Nitrogen (mg/dL)", min_value=0.0, value=None, step=0.1)
            uric_acid = st.number_input("Serum Uric Acid (mg/dL)", min_value=0.0, value=None, step=0.1)
        
        # Lipid Panel
        with st.expander("🧴 Lipid Panel", expanded=True):
            cholesterol = st.number_input("Total Cholesterol (mg/dL)", min_value=0.0, value=None, step=0.1)
            hdl = st.number_input("HDL-C (mg/dL)", min_value=0.0, value=None, step=0.1)
            ldl = st.number_input("LDL-C (mg/dL)", min_value=0.0, value=None, step=0.1)
            triglycerides = st.number_input("Triglycerides (mg/dL)", min_value=0.0, value=None, step=0.1)
            vldl = st.number_input("VLDL-C (mg/dL)", min_value=0.0, value=None, step=0.1)
        
        # Liver Function
        with st.expander("🫀 Liver Function Tests"):
            alt = st.number_input("ALT (U/L)", min_value=0.0, value=None, step=0.1)
            ast = st.number_input("AST (U/L)", min_value=0.0, value=None, step=0.1)
            alp = st.number_input("ALP (IU/L)", min_value=0.0, value=None, step=0.1)
        
        # Additional CKD markers
        with st.expander("🫘 Kidney Function"):
            calcium = st.number_input("Serum Calcium (mg/dL)", min_value=0.0, value=None, step=0.1)
            phosphate = st.number_input("Serum Phosphate (mg/dL)", min_value=0.0, value=None, step=0.1)
        
        if st.button("🔬 Perform AIKEGG Analysis", type="primary"):
            lab_values = {}
            
            # Collect all non-None values
            values_map = {
                'glucose': glucose, 'bun': bun, 'creatinine': creatinine,
                'cholesterol': cholesterol, 'hdl': hdl, 'ldl': ldl,
                'triglycerides': triglycerides, 'vldl': vldl,
                'alt': alt, 'ast': ast, 'alp': alp, 'uric_acid': uric_acid,
                'calcium': calcium, 'phosphate': phosphate
            }
            
            for key, value in values_map.items():
                if value is not None:
                    lab_values[key] = value
            
            if not lab_values:
                st.error("Please enter at least one lab value to analyze.")
                return
            
            pathway_scores = get_pathway_perturbation_score(lab_values, metabolic_db, ref_ranges)
            sorted_pathways = sorted(pathway_scores.items(), key=lambda x: x[1]['total_score'], reverse=True)
            
            # Get affected pathways analysis
            directly_affected, at_risk = get_affected_pathways_analysis(lab_values, metabolic_db, ref_ranges)
            
            st.session_state['analysis_results'] = {
                'lab_values': lab_values,
                'sorted_pathways': sorted_pathways,
                'directly_affected': directly_affected,
                'at_risk': at_risk,
                'patient_info': {
                    'patient_id': patient_id,
                    'age': age,
                    'gender': gender
                }
            }
    
    with col2:
        if 'analysis_results' in st.session_state:
            results = st.session_state['analysis_results']
            lab_values = results['lab_values']
            sorted_pathways = results['sorted_pathways']
            directly_affected = results['directly_affected']
            at_risk = results['at_risk']
            patient_info = results['patient_info']
            
            tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["Overview", "Affected Pathways", "Pathway Details", "Network View", "Thermodynamics & Flux", "Enhanced Report"])
            
            with tab1:
                st.subheader("🎯 Analysis Overview")
                
                # Key metrics
                col_a, col_b, col_c, col_d = st.columns(4)
                with col_a:
                    st.metric("Parameters", len(lab_values))
                with col_b:
                    st.metric("Directly Affected", len(directly_affected), delta=f"{len(directly_affected)} pathways")
                with col_c:
                    st.metric("At Risk", len(at_risk), delta=f"{len(at_risk)} pathways")
                with col_d:
                    st.metric("Top Score", f"{sorted_pathways[0][1]['total_score']:.2f}")
                
                # Overview chart
                st.plotly_chart(create_overview_chart(sorted_pathways), use_container_width=True)
                
                # Quick summary
                st.subheader("🚨 Immediate Attention Required")
                if directly_affected:
                    for pathway in directly_affected[:3]:
                        st.error(f"**{pathway['name']}** ({pathway['disease']}) - Score: {pathway['score']:.2f}")
                else:
                    st.success("No pathways require immediate attention")
            
            with tab2:
                st.subheader("📊 Pathway Perturbation Analysis")
               # Create Excel-like table
     
                table_data = create_pathway_table(sorted_pathways, directly_affected, at_risk)
                
                df = pd.DataFrame(table_data)

                def highlight_status(val):
                    if '🔴' in val:
                        return 'background-color: #ffebee; color: #c62828; font-weight: bold'
                    elif '🟡' in val:
                        return 'background-color: #fff3e0; color: #ef6c00; font-weight: bold'
                    elif '🟢' in val:
                        return 'background-color: #e8f5e8; color: #2e7d32; font-weight: bold'
                    return ''
                def highlight_scores(val):
                    try:
                        score = float(val)
                        if score > 2.0:
                            return 'background-color: #ffcdd2; font-weight: bold'
                        elif score > 1.0:
                            return 'background-color: #ffe0b2; font-weight: bold'
                        return ''
                    except:
                        return ''
                styled_df = df.style.applymap(highlight_status, subset=['Status']) \
                                        .applymap(highlight_scores, subset=['Total Score']) \
                                        .set_properties(**{
                                            'font-size': '12px',
                                            'border': '1px solid #ddd'
                                            })
                st.dataframe(styled_df, use_container_width=True, height=400)
                # Summary statistics
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Total Pathways", len(sorted_pathways))
                with col2:
                    st.metric("Directly Affected", len(directly_affected))
                with col3:
                    st.metric("At Risk", len(at_risk))
                with col4:
                    avg_score = sum(data[1]['total_score'] for data in sorted_pathways) / len(sorted_pathways)
                    st.metric("Average Score", f"{avg_score:.2f}")    

                if directly_affected:
                    st.error("🔴 **DIRECTLY AFFECTED PATHWAYS**")
                    for pathway in directly_affected:
                        with st.expander(f"{pathway['name']} ({pathway['disease']}) - Score: {pathway['score']:.2f}"):
                            st.write("**Abnormal Biomarkers:**")
                            if pathway['affected_biomarkers']:
                                for biomarker in pathway['affected_biomarkers']:
                                    st.write(f"• {biomarker}")
                            else:
                                st.write("• Based on clinical correlation and pathway relevance")

                if at_risk:
                    st.warning("🟡 **PATHWAYS AT RISK**")
                    for pathway in at_risk:
                        with st.expander(f"{pathway['name']} ({pathway['disease']}) - Score: {pathway['score']:.2f}"):
                            st.write("**Risk Factors:**")
                            if pathway['risk_factors']:
                                for factor in pathway['risk_factors']:
                                    st.write(f"• {factor}")
                            else:
                                st.write("• Moderate perturbation detected")                                        
            
            with tab3:
                st.subheader("🔬 Detailed Pathway Information")
                
                # Select pathway for detailed view
                pathway_options = {f"{data[1]['pathway']['name']} ({data[1]['pathway']['disease']})": data[0] 
                                 for data in sorted_pathways}
                selected_pathway_name = st.selectbox("Select pathway for detailed analysis:", 
                                                   list(pathway_options.keys()))
                
                if selected_pathway_name:
                    selected_pathway_id = pathway_options[selected_pathway_name]
                    pathway_data = {pid: data for pid, data in sorted_pathways}
                    display_pathway_details(pathway_data, selected_pathway_id)
            
            with tab4:
                st.subheader("🕸️ KEGG-Style Metabolic Pathway Network")
                
                # Add pathway selector for detailed view
                col_select, col_button = st.columns([3, 1])
                
                with col_select:
                    pathway_options = ["Overview"] + [f"{data[1]['pathway']['name']} ({data[1]['pathway']['disease']})" 
                                                for data in sorted_pathways]
                    selected_view = st.selectbox("Select view:", pathway_options)
                
                with col_button:
                    st.write("")  # Spacing
                    show_overview = st.button("🔄 Show Overview")
                
                pathway_data = {pid: data for pid, data in sorted_pathways}
                
                if show_overview or selected_view == "Overview":
                    network_fig = create_kegg_style_pathway_chart(pathway_data)
                    st.plotly_chart(network_fig, use_container_width=True)
                    st.info("🔵 Blue circles = Compounds | 🟩 Green squares = Enzymes | ➡️ Red arrows = Reactions")
                else:
                    # Find selected pathway ID
                    selected_pathway_id = None
                    for pid, data in sorted_pathways:
                        pathway_name = f"{data['pathway']['name']} ({data['pathway']['disease']})"
                        if pathway_name == selected_view:
                            selected_pathway_id = pid
                            break
                    
                    if selected_pathway_id:
                        detailed_fig = create_kegg_style_pathway_chart(pathway_data, selected_pathway_id)
                        st.plotly_chart(detailed_fig, use_container_width=True)
                        
                        # Show additional pathway details
                        pathway_info = pathway_data[selected_pathway_id]['pathway']
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.write("**Key Metabolites:**")
                            for comp_name in pathway_info['compound_names'][:4]:
                                st.write(f"• {comp_name}")
                        
                        with col2:
                            st.write("**Key Enzymes:**")
                            for enz_name in pathway_info['enzyme_names'][:3]:
                                st.write(f"• {enz_name}")
                        
                        st.write("**Biochemical Reactions:**")
                        for reaction_name in pathway_info['reaction_names'][:3]:
                            st.write(f"➡️ {reaction_name}")
            
            with tab5:
                
                st.info("This analysis shows the thermodynamic favorability and metabolic flux perturbations")
                
                # Thermodynamic Analysis
                st.subheader("🔥 Thermodynamic Favorability")
                col_thermo, col_info = st.columns([2, 1])
                
                with col_thermo:
                    st.plotly_chart(create_thermodynamics_chart(sorted_pathways), use_container_width=True)
                
                with col_info:
                    st.write("**Interpretation:**")
                    st.write("• ΔG < -50: Highly Favorable")
                    st.write("• ΔG -50 to -20: Moderately Favorable")
                    st.write("• ΔG > -20: Less Favorable")
                    
                    # Show thermodynamic summary table
                    thermo_data = []
                    for pathway_id, data in sorted_pathways[:5]:
                        favorability = "Highly Favorable" if data['pathway']['deltaG'] < -50 else "Moderately Favorable" if data['pathway']['deltaG'] < -20 else "Less Favorable"
                        thermo_data.append({
                            'Pathway': data['pathway']['name'][:20] + '...' if len(data['pathway']['name']) > 20 else data['pathway']['name'],
                            'ΔG (kJ/mol)': data['pathway']['deltaG'],
                            'Favorability': favorability
                        })
                    
                    st.write("**Top 5 Pathways:**")
                    st.dataframe(pd.DataFrame(thermo_data), use_container_width=True)
                
                st.divider()
                
                # Flux Analysis
                st.subheader("📊 Metabolic Flux Perturbation")
                col_flux, col_flux_info = st.columns([2, 1])
                
                with col_flux:
                    st.plotly_chart(create_flux_chart(sorted_pathways), use_container_width=True)
                
                with col_flux_info:
                    st.write("**Flux Analysis:**")
                    st.write("• Red area: Perturbed flux state")
                    st.write("• Green area: Normal flux state")
                    st.write("• Larger deviation = Higher perturbation")
                    
                    # Flux perturbation summary
                    flux_data = []
                    for pathway_id, data in sorted_pathways[:5]:
                        perturbation_percent = data['perturbation_score'] * 30  # Convert to percentage
                        flux_data.append({
                            'Pathway': data['pathway']['name'][:20] + '...' if len(data['pathway']['name']) > 20 else data['pathway']['name'],
                            'Perturbation %': f"{perturbation_percent:.1f}%",
                            'Status': "High" if perturbation_percent > 45 else "Moderate" if perturbation_percent > 15 else "Low"
                        })
                    
                    st.write("**Flux Perturbations:**")
                    st.dataframe(pd.DataFrame(flux_data), use_container_width=True)
            
            with tab6:
                st.subheader("📋 Enhanced Clinical Report")
                
                enhanced_report = generate_enhanced_report_content(lab_values, sorted_pathways, patient_info, ref_ranges)
                
                st.text_area("Clinical Report", enhanced_report, height=600)
                
                col_download1, col_download2, col_actions = st.columns(3)
                
                with col_download1:
                    st.download_button(
                        label="📄 Download TXT Report",
                        data=enhanced_report,
                        file_name=f"AIKEGG_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                        mime="text/plain"
                    )
                
                with col_download2:
                    try:
                        pdf_buffer = generate_pdf_report(lab_values, sorted_pathways, patient_info, ref_ranges)
                        st.download_button(
                            label="📑 Download PDF Report",
                            data=pdf_buffer,
                            file_name=f"AIKEGG_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                            mime="application/pdf"
                        )
                    except Exception as e:
                        st.error(f"PDF generation failed. Install reportlab: pip install reportlab")
                
                

if __name__ == "__main__":
    main()