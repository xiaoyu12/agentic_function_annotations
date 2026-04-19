$ErrorActionPreference = 'Stop'

$paperDir = Split-Path -Parent $PSScriptRoot
$projectDir = Split-Path -Parent $paperDir
$comparisonDir = Join-Path $projectDir 'comparison_merged_annotation'
$figDir = Join-Path $paperDir 'figures'
$suppDir = Join-Path $paperDir 'supplement'
New-Item -ItemType Directory -Force $figDir | Out-Null
New-Item -ItemType Directory -Force $suppDir | Out-Null

Add-Type -AssemblyName System.Drawing

function New-Color {
    param([string]$Hex)
    return [System.Drawing.ColorTranslator]::FromHtml($Hex)
}

function New-SolidBrush {
    param([string]$Hex)
    return [System.Drawing.SolidBrush]::new((New-Color $Hex))
}

function Draw-Text {
    param(
        [System.Drawing.Graphics]$Graphics,
        [string]$Text,
        [System.Drawing.Font]$Font,
        [System.Drawing.Brush]$Brush,
        [float]$X,
        [float]$Y
    )
    $Graphics.DrawString($Text, $Font, $Brush, $X, $Y)
}

function Draw-CenteredText {
    param(
        [System.Drawing.Graphics]$Graphics,
        [string]$Text,
        [System.Drawing.Font]$Font,
        [System.Drawing.Brush]$Brush,
        [System.Drawing.RectangleF]$Rectangle
    )
    $format = [System.Drawing.StringFormat]::new()
    $format.Alignment = [System.Drawing.StringAlignment]::Center
    $format.LineAlignment = [System.Drawing.StringAlignment]::Center
    $Graphics.DrawString($Text, $Font, $Brush, $Rectangle, $format)
    $format.Dispose()
}

function Save-Bitmap {
    param(
        [System.Drawing.Bitmap]$Bitmap,
        [string]$Path
    )
    $Bitmap.Save($Path, [System.Drawing.Imaging.ImageFormat]::Png)
}

function Parse-LabelCounts {
    param([string]$Labels)
    $counts = @{
        High = 0
        Moderate = 0
        Watchlist = 0
        Low = 0
    }
    foreach ($part in ($Labels -split ';')) {
        $trimmed = $part.Trim()
        if ($trimmed -notmatch '=') { continue }
        $pieces = $trimmed -split '=', 2
        $key = $pieces[0].Trim().ToLowerInvariant()
        $value = [int]$pieces[1].Trim()
        if ($key.StartsWith('high')) { $counts.High = $value }
        elseif ($key.StartsWith('moderate')) { $counts.Moderate = $value }
        elseif ($key.StartsWith('watchlist')) { $counts.Watchlist = $value }
        elseif ($key.StartsWith('low')) { $counts.Low = $value }
    }
    return $counts
}

function Get-RunDisplay {
    param([string]$Run)
    if ($Run -match '^Claude_app_run(\d)$') { return "Claude App R$($Matches[1])" }
    if ($Run -match '^Claude_code_run(\d)$') { return "Claude Code R$($Matches[1])" }
    if ($Run -match '^Codex_app_run(\d)$') { return "Codex App R$($Matches[1])" }
    return $Run
}

function Escape-Latex {
    param([AllowNull()][string]$Text)
    if ([string]::IsNullOrEmpty($Text)) { return '' }
    $s = $Text
    $s = $s.Replace('\', '\textbackslash{}')
    $s = $s.Replace('&', '\&')
    $s = $s.Replace('%', '\%')
    $s = $s.Replace('$', '\$')
    $s = $s.Replace('#', '\#')
    $s = $s.Replace('_', '\_')
    $s = $s.Replace('{', '\{')
    $s = $s.Replace('}', '\}')
    $s = $s.Replace('^', '\textasciicircum{}')
    $s = $s.Replace('~', '\textasciitilde{}')
    return $s
}

function Clamp01 {
    param([double]$Value)
    if ($Value -lt 0) { return 0.0 }
    if ($Value -gt 1) { return 1.0 }
    return $Value
}

function Interpolate-Blue {
    param([double]$Value)
    $v = Clamp01 $Value
    $r1 = 247; $g1 = 251; $b1 = 255
    $r2 = 8; $g2 = 48; $b2 = 107
    $r = [int]($r1 + ($r2 - $r1) * $v)
    $g = [int]($g1 + ($g2 - $g1) * $v)
    $b = [int]($b1 + ($b2 - $b1) * $v)
    return [System.Drawing.Color]::FromArgb($r, $g, $b)
}

$methodPath = Join-Path $comparisonDir 'method_consistency_summary.tsv'
$finalPath = Join-Path $comparisonDir 'final_orthogroup_annotations.tsv'
$methodRows = Import-Csv -Delimiter "`t" $methodPath
$finalRows = Import-Csv -Delimiter "`t" $finalPath

# Figure 1: final curated relevance distribution.
$finalOrder = @('High', 'Moderate', 'Watchlist', 'Low')
$finalPalette = @{
    Low = '#999999'
    Watchlist = '#56B4E9'
    Moderate = '#E69F00'
    High = '#D55E00'
}
$finalCounts = @{}
foreach ($tier in $finalOrder) {
    $finalCounts[$tier] = @($finalRows | Where-Object { $_.final_calcification_relevance -eq $tier }).Count
}
$width = 1600
$height = 900
$bmp = [System.Drawing.Bitmap]::new($width, $height)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.Clear([System.Drawing.Color]::White)
$black = New-SolidBrush '#111111'
$gray = New-SolidBrush '#666666'
$titleFont = [System.Drawing.Font]::new('Arial', 34, [System.Drawing.FontStyle]::Bold)
$axisFont = [System.Drawing.Font]::new('Arial', 24)
$labelFont = [System.Drawing.Font]::new('Arial', 22)
$smallFont = [System.Drawing.Font]::new('Arial', 18)
$whiteBrush = New-SolidBrush '#FFFFFF'
Draw-Text $g 'Final merged calcification-relevance tiers' $titleFont $black 90 65
Draw-Text $g 'Curated consensus after evidence audit across 73 orthogroups' $smallFont $gray 90 112
$plotLeft = 360
$plotRight = 1340
$plotTop = 220
$barHeight = 90
$gap = 50
$axisY = $plotTop + $finalOrder.Count * ($barHeight + $gap) - 10
$maxCount = 45.0
$plotWidth = $plotRight - $plotLeft
$axisPen = [System.Drawing.Pen]::new((New-Color '#222222'), 2)
$gridPen = [System.Drawing.Pen]::new((New-Color '#DDDDDD'), 1)
foreach ($tick in @(0, 10, 20, 30, 40)) {
    $x = $plotLeft + ($tick / $maxCount) * $plotWidth
    $g.DrawLine($gridPen, [float]$x, [float]($plotTop - 25), [float]$x, [float]($axisY - 15))
    $g.DrawLine($axisPen, [float]$x, [float]($axisY - 15), [float]$x, [float]$axisY)
    Draw-CenteredText $g ([string]$tick) $smallFont $black ([System.Drawing.RectangleF]::new([float]($x - 35), [float]($axisY + 10), 70, 32))
}
$g.DrawLine($axisPen, $plotLeft, [float]$axisY, $plotRight, [float]$axisY)
for ($i = 0; $i -lt $finalOrder.Count; $i++) {
    $tier = $finalOrder[$i]
    $count = [int]$finalCounts[$tier]
    $pct = 100.0 * $count / 73.0
    $y = $plotTop + $i * ($barHeight + $gap)
    Draw-Text $g $tier $axisFont $black 115 ($y + 28)
    $barW = [float](($count / $maxCount) * $plotWidth)
    $brush = New-SolidBrush $finalPalette[$tier]
    $g.FillRectangle($brush, [float]$plotLeft, [float]$y, $barW, [float]$barHeight)
    $g.DrawRectangle([System.Drawing.Pens]::Black, $plotLeft, $y, [int]$barW, $barHeight)
    $label = [string]::Format('{0} ({1:0.0}%)', $count, $pct)
    if ($barW -gt 155) {
        $textBrush = if ($tier -in @('High', 'Moderate')) { $whiteBrush } else { $black }
        Draw-CenteredText $g $label $labelFont $textBrush ([System.Drawing.RectangleF]::new([float]$plotLeft, [float]$y, $barW, [float]$barHeight))
    } else {
        Draw-Text $g $label $labelFont $black ([float]($plotLeft + $barW + 18)) ($y + 29)
    }
    $brush.Dispose()
}
Draw-CenteredText $g 'Number of orthogroups' $axisFont $black ([System.Drawing.RectangleF]::new([float]$plotLeft, [float]($axisY + 52), [float]$plotWidth, 45))
Save-Bitmap $bmp (Join-Path $figDir 'final_relevance_distribution.png')
$g.Dispose()
$bmp.Dispose()

# Figure 2: stacked relevance distributions by run.
$runData = foreach ($row in $methodRows) {
    $counts = Parse-LabelCounts $row.labels
    [pscustomobject]@{
        Run = $row.run
        Display = Get-RunDisplay $row.run
        High = $counts.High
        Moderate = $counts.Moderate
        Watchlist = $counts.Watchlist
        Low = $counts.Low
        MeanScore = [double]$row.mean_score
        ScoreSd = [double]$row.score_sd
        Missing = [int]$row.missing_ogs
    }
}

$width = 1900
$height = 1220
$bmp = [System.Drawing.Bitmap]::new($width, $height)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.Clear([System.Drawing.Color]::White)
$black = New-SolidBrush '#111111'
$gray = New-SolidBrush '#666666'
$titleFont = [System.Drawing.Font]::new('Arial', 34, [System.Drawing.FontStyle]::Bold)
$axisFont = [System.Drawing.Font]::new('Arial', 22)
$smallFont = [System.Drawing.Font]::new('Arial', 18)
$whiteBrush = New-SolidBrush '#FFFFFF'
$palette = @{
    Low = '#999999'
    Watchlist = '#56B4E9'
    Moderate = '#E69F00'
    High = '#D55E00'
}
Draw-Text $g 'Normalized calcification-relevance calls by agent run' $titleFont $black 110 55
Draw-Text $g 'All bars sum to 73 orthogroups; colors show normalized evidence tier.' $smallFont $gray 110 100
$left = 350
$right = 110
$top = 230
$barHeight = 58
$gap = 42
$plotWidth = $width - $left - $right
$scale = $plotWidth / 73.0
$order = @('Low', 'Watchlist', 'Moderate', 'High')
for ($i = 0; $i -lt $runData.Count; $i++) {
    $row = $runData[$i]
    $y = $top + $i * ($barHeight + $gap)
    Draw-Text $g $row.Display $axisFont $black 70 ($y + 14)
    $x = $left
    foreach ($tier in $order) {
        $count = [int]$row.$tier
        $segmentWidth = [float]($count * $scale)
        $brush = New-SolidBrush $palette[$tier]
        $g.FillRectangle($brush, [float]$x, [float]$y, $segmentWidth, [float]$barHeight)
        $g.DrawRectangle([System.Drawing.Pens]::White, [int]$x, [int]$y, [int]$segmentWidth, [int]$barHeight)
        if ($segmentWidth -gt 42 -and $count -gt 0) {
            $rect = [System.Drawing.RectangleF]::new([float]$x, [float]$y, $segmentWidth, [float]$barHeight)
            $textBrush = if ($tier -eq 'Low' -or $tier -eq 'Watchlist') { $black } else { $whiteBrush }
            Draw-CenteredText $g ([string]$count) $smallFont $textBrush $rect
        }
        $brush.Dispose()
        $x += $segmentWidth
    }
}

$axisY = $top + $runData.Count * ($barHeight + $gap) - 18
$pen = [System.Drawing.Pen]::new((New-Color '#333333'), 2)
$g.DrawLine($pen, $left, $axisY, $left + $plotWidth, $axisY)
foreach ($tick in @(0, 20, 40, 60, 73)) {
    $x = $left + $tick * $scale
    $g.DrawLine($pen, [float]$x, [float]$axisY, [float]$x, [float]($axisY + 12))
    Draw-CenteredText $g ([string]$tick) $smallFont $black ([System.Drawing.RectangleF]::new([float]($x - 35), [float]($axisY + 15), 70, 35))
}
Draw-CenteredText $g 'Number of orthogroups' $axisFont $black ([System.Drawing.RectangleF]::new([float]$left, [float]($axisY + 55), [float]$plotWidth, 40))

$legendX = 1060
$legendY = 132
foreach ($tier in @('High', 'Moderate', 'Watchlist', 'Low')) {
    $brush = New-SolidBrush $palette[$tier]
    $g.FillRectangle($brush, $legendX, $legendY, 32, 24)
    $g.DrawRectangle([System.Drawing.Pens]::Black, $legendX, $legendY, 32, 24)
    Draw-Text $g $tier $smallFont $black ($legendX + 45) ($legendY - 3)
    $legendX += 205
    $brush.Dispose()
}
Save-Bitmap $bmp (Join-Path $figDir 'run_relevance_distribution.png')
$g.Dispose()
$bmp.Dispose()

# Figure 3: within-method exact agreement dot plot.
$agreement = @(
    [pscustomobject]@{ Method = 'Claude App'; Pair = 'R1-R2'; Value = 0.397 },
    [pscustomobject]@{ Method = 'Claude App'; Pair = 'R1-R3'; Value = 0.603 },
    [pscustomobject]@{ Method = 'Claude App'; Pair = 'R2-R3'; Value = 0.685 },
    [pscustomobject]@{ Method = 'Claude Code + skills'; Pair = 'R1-R2'; Value = 0.466 },
    [pscustomobject]@{ Method = 'Claude Code + skills'; Pair = 'R1-R3'; Value = 0.342 },
    [pscustomobject]@{ Method = 'Claude Code + skills'; Pair = 'R2-R3'; Value = 0.740 },
    [pscustomobject]@{ Method = 'Codex App + skills'; Pair = 'R1-R2'; Value = 0.630 },
    [pscustomobject]@{ Method = 'Codex App + skills'; Pair = 'R1-R3'; Value = 0.411 },
    [pscustomobject]@{ Method = 'Codex App + skills'; Pair = 'R2-R3'; Value = 0.575 }
)
$width = 1500
$height = 1000
$bmp = [System.Drawing.Bitmap]::new($width, $height)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.Clear([System.Drawing.Color]::White)
$titleFont = [System.Drawing.Font]::new('Arial', 34, [System.Drawing.FontStyle]::Bold)
$axisFont = [System.Drawing.Font]::new('Arial', 22)
$smallFont = [System.Drawing.Font]::new('Arial', 18)
Draw-Text $g 'Run-to-run consistency within each method' $titleFont $black 90 60
Draw-Text $g 'Exact agreement of normalized relevance tiers across 73 orthogroups.' $smallFont $gray 90 108
$plotLeft = 170
$plotTop = 170
$plotRight = 1320
$plotBottom = 800
$yMin = 0.30
$yMax = 0.80
$gridPen = [System.Drawing.Pen]::new((New-Color '#DDDDDD'), 1)
$axisPen = [System.Drawing.Pen]::new((New-Color '#222222'), 2)
for ($tick = 0.30; $tick -le 0.801; $tick += 0.10) {
    $y = $plotBottom - (($tick - $yMin) / ($yMax - $yMin)) * ($plotBottom - $plotTop)
    $g.DrawLine($gridPen, $plotLeft, [float]$y, $plotRight, [float]$y)
    Draw-Text $g ([string]::Format('{0:0.0}', $tick)) $smallFont $black 88 ([float]($y - 14))
}
$g.DrawLine($axisPen, $plotLeft, $plotTop, $plotLeft, $plotBottom)
$g.DrawLine($axisPen, $plotLeft, $plotBottom, $plotRight, $plotBottom)
$methodX = @{
    'Claude App' = 360
    'Claude Code + skills' = 750
    'Codex App + skills' = 1140
}
$methodColors = @{
    'Claude App' = '#0072B2'
    'Claude Code + skills' = '#009E73'
    'Codex App + skills' = '#D55E00'
}
foreach ($method in $methodX.Keys) {
    $x = $methodX[$method]
    Draw-CenteredText $g $method $axisFont $black ([System.Drawing.RectangleF]::new([float]($x - 170), [float]($plotBottom + 35), 340, 60))
    $vals = @($agreement | Where-Object { $_.Method -eq $method } | ForEach-Object { $_.Value })
    $mean = ($vals | Measure-Object -Average).Average
    $meanY = $plotBottom - (($mean - $yMin) / ($yMax - $yMin)) * ($plotBottom - $plotTop)
    $meanPen = [System.Drawing.Pen]::new((New-Color $methodColors[$method]), 5)
    $g.DrawLine($meanPen, [float]($x - 85), [float]$meanY, [float]($x + 85), [float]$meanY)
    $meanPen.Dispose()
}
$jitter = @{ 'R1-R2' = -65; 'R1-R3' = 0; 'R2-R3' = 65 }
foreach ($point in $agreement) {
    $x = $methodX[$point.Method] + $jitter[$point.Pair]
    $y = $plotBottom - (($point.Value - $yMin) / ($yMax - $yMin)) * ($plotBottom - $plotTop)
    $brush = New-SolidBrush $methodColors[$point.Method]
    $g.FillEllipse($brush, [float]($x - 18), [float]($y - 18), 36, 36)
    $g.DrawEllipse([System.Drawing.Pens]::Black, [float]($x - 18), [float]($y - 18), 36, 36)
    Draw-CenteredText $g $point.Pair $smallFont $black ([System.Drawing.RectangleF]::new([float]($x - 45), [float]($y - 55), 90, 25))
    Draw-CenteredText $g ([string]::Format('{0:0.000}', $point.Value)) $smallFont $black ([System.Drawing.RectangleF]::new([float]($x - 55), [float]($y + 22), 110, 28))
    $brush.Dispose()
}
$state = $g.Save()
$g.TranslateTransform(42, [float](($plotTop + $plotBottom) / 2 + 145))
$g.RotateTransform(-90)
Draw-CenteredText $g 'Agreement fraction' $axisFont $black ([System.Drawing.RectangleF]::new(0, 0, 360, 40))
$g.Restore($state)
Save-Bitmap $bmp (Join-Path $figDir 'within_method_agreement.png')
$g.Dispose()
$bmp.Dispose()

# Figure 4: evidence heatmap for final high and moderate candidates.
$candidates = $finalRows |
    Where-Object { $_.final_calcification_relevance -in @('High', 'Moderate') } |
    Sort-Object @{ Expression = { if ($_.final_calcification_relevance -eq 'High') { 0 } else { 1 } } }, orthogroup
$features = @(
    @{ Key = 'signalp'; Label = 'SignalP'; Scale = 1.0 },
    @{ Key = 'tm'; Label = 'TM'; Scale = 1.0 },
    @{ Key = 'blast'; Label = 'BLAST'; Scale = 1.0 },
    @{ Key = 'pfam'; Label = 'Pfam'; Scale = 1.0 },
    @{ Key = 'mean'; Label = 'Model'; Scale = 4.0 },
    @{ Key = 'sd'; Label = 'Disagr.'; Scale = 1.5 }
)
$matrix = foreach ($row in $candidates) {
    $n = [double]$row.sequence_count
    [pscustomobject]@{
        Orthogroup = $row.orthogroup
        Tier = $row.final_calcification_relevance
        signalp = ([double]$row.signalp_positive / $n)
        tm = ([double]$row.deeptmhmm_tm_or_beta_count / $n)
        blast = ([double]$row.blast_hit_queries / $n)
        pfam = ([double]$row.pfam_hit_queries / $n)
        mean = ([double]$row.model_consensus_mean_score / 4.0)
        sd = ([double]$row.model_disagreement_sd / 1.5)
        signalpRaw = ([double]$row.signalp_positive / $n)
        tmRaw = ([double]$row.deeptmhmm_tm_or_beta_count / $n)
        blastRaw = ([double]$row.blast_hit_queries / $n)
        pfamRaw = ([double]$row.pfam_hit_queries / $n)
        meanRaw = [double]$row.model_consensus_mean_score
        sdRaw = [double]$row.model_disagreement_sd
    }
}

$width = 1750
$height = 1280
$bmp = [System.Drawing.Bitmap]::new($width, $height)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.Clear([System.Drawing.Color]::White)
$titleFont = [System.Drawing.Font]::new('Arial', 32, [System.Drawing.FontStyle]::Bold)
$axisFont = [System.Drawing.Font]::new('Arial', 22)
$smallFont = [System.Drawing.Font]::new('Arial', 17)
$tinyFont = [System.Drawing.Font]::new('Arial', 15)
Draw-Text $g 'Evidence profiles of final high and moderate candidates' $titleFont $black 80 55
Draw-Text $g 'Fractions are sequence-level coverage except Model score and Disagreement, which are scaled for color only.' $smallFont $gray 80 102
$left = 370
$top = 205
$cellW = 170
$cellH = 68
for ($c = 0; $c -lt $features.Count; $c++) {
    $x = $left + $c * $cellW
    Draw-CenteredText $g $features[$c].Label $axisFont $black ([System.Drawing.RectangleF]::new([float]$x, 150, [float]$cellW, 45))
}
for ($r = 0; $r -lt $matrix.Count; $r++) {
    $row = $matrix[$r]
    $y = $top + $r * $cellH
    $tierMarker = if ($row.Tier -eq 'High') { 'H' } else { 'M' }
    Draw-Text $g ("$($row.Orthogroup) ($tierMarker)") $axisFont $black 70 ($y + 17)
    for ($c = 0; $c -lt $features.Count; $c++) {
        $feature = $features[$c]
        $x = $left + $c * $cellW
        $scaled = Clamp01 ([double]$row.($feature.Key))
        $brush = [System.Drawing.SolidBrush]::new((Interpolate-Blue $scaled))
        $g.FillRectangle($brush, [float]$x, [float]$y, [float]($cellW - 4), [float]($cellH - 4))
        $g.DrawRectangle([System.Drawing.Pens]::White, $x, $y, $cellW - 4, $cellH - 4)
        $rawKey = "$($feature.Key)Raw"
        $raw = [double]$row.$rawKey
        $text = if ($feature.Key -in @('mean', 'sd')) { [string]::Format('{0:0.00}', $raw) } else { [string]::Format('{0:0.00}', $raw) }
        $textBrush = if ($scaled -gt 0.55) { $whiteBrush } else { $black }
        Draw-CenteredText $g $text $smallFont $textBrush ([System.Drawing.RectangleF]::new([float]$x, [float]$y, [float]($cellW - 4), [float]($cellH - 4)))
        $brush.Dispose()
    }
}
$legendLeft = $left
$legendTop = $top + $matrix.Count * $cellH + 50
for ($i = 0; $i -le 100; $i++) {
    $color = Interpolate-Blue ($i / 100.0)
    $penColor = [System.Drawing.Pen]::new($color, 3)
    $x = $legendLeft + $i * 5
    $g.DrawLine($penColor, $x, $legendTop, $x, $legendTop + 35)
    $penColor.Dispose()
}
$g.DrawRectangle([System.Drawing.Pens]::Black, $legendLeft, $legendTop, 500, 35)
Draw-Text $g 'Low' $smallFont $black ($legendLeft - 5) ($legendTop + 42)
Draw-Text $g 'High' $smallFont $black ($legendLeft + 455) ($legendTop + 42)
Draw-Text $g 'H = final high; M = final moderate' $smallFont $gray ($legendLeft + 640) ($legendTop + 10)
Save-Bitmap $bmp (Join-Path $figDir 'candidate_evidence_heatmap.png')
$g.Dispose()
$bmp.Dispose()

# Supplementary Table S1: complete final annotation table.
$finalTablePath = Join-Path $suppDir 'final_annotations_table.tex'
$lines = [System.Collections.Generic.List[string]]::new()
$lines.Add('\begingroup')
$lines.Add('\footnotesize')
$lines.Add('\setlength{\tabcolsep}{3pt}')
$lines.Add('\begin{longtable}{p{1.9cm}r p{4.55cm}p{1.8cm}p{5.55cm}}')
$lines.Add('\caption{Complete final merged annotations for the 73 orthogroups. The table is derived from \texttt{final\_orthogroup\_annotations.tsv}.}\label{tab:s_final_annotations}\\')
$lines.Add('\toprule')
$lines.Add('Orthogroup & n & Final function & Relevance & Curated rationale\\')
$lines.Add('\midrule')
$lines.Add('\endfirsthead')
$lines.Add('\caption[]{Complete final merged annotations for the 73 orthogroups (continued).}\\')
$lines.Add('\toprule')
$lines.Add('Orthogroup & n & Final function & Relevance & Curated rationale\\')
$lines.Add('\midrule')
$lines.Add('\endhead')
$lines.Add('\midrule')
$lines.Add('\multicolumn{5}{r}{Continued on next page}\\')
$lines.Add('\endfoot')
$lines.Add('\bottomrule')
$lines.Add('\endlastfoot')
foreach ($row in $finalRows) {
    $og = Escape-Latex $row.orthogroup
    $n = Escape-Latex $row.sequence_count
    $func = Escape-Latex $row.final_function
    $rel = Escape-Latex $row.final_calcification_relevance
    $rat = Escape-Latex $row.final_rationale
    $lines.Add("$og & $n & $func & $rel & $rat\\")
}
$lines.Add('\end{longtable}')
$lines.Add('\endgroup')
Set-Content -Path $finalTablePath -Value $lines -Encoding UTF8

# Supplementary Table S2: method summary.
$methodTablePath = Join-Path $suppDir 'method_summary_table.tex'
$lines = [System.Collections.Generic.List[string]]::new()
$lines.Add('\begin{table}[ht]')
$lines.Add('\centering')
$lines.Add('\small')
$lines.Add('\caption{Full method-level consistency summary. Label distributions were normalized from each agent output.}')
$lines.Add('\label{tab:s_method_summary}')
$lines.Add('\begin{tabularx}{\linewidth}{lYrrr}')
$lines.Add('\toprule')
$lines.Add('Run & Normalized label distribution & Mean score & Score SD & Missing OGs\\')
$lines.Add('\midrule')
foreach ($row in $methodRows) {
    $run = Escape-Latex $row.run
    $labels = Escape-Latex $row.labels
    $mean = Escape-Latex $row.mean_score
    $sd = Escape-Latex $row.score_sd
    $missing = Escape-Latex $row.missing_ogs
    $lines.Add("\texttt{$run} & $labels & $mean & $sd & $missing\\")
}
$lines.Add('\bottomrule')
$lines.Add('\end{tabularx}')
$lines.Add('\end{table}')
Set-Content -Path $methodTablePath -Value $lines -Encoding UTF8

# Supplementary Table S3: pairwise agreement values.
$pairwiseTablePath = Join-Path $suppDir 'pairwise_agreement_table.tex'
$lines = [System.Collections.Generic.List[string]]::new()
$lines.Add('\begin{table}[ht]')
$lines.Add('\centering')
$lines.Add('\small')
$lines.Add('\caption{Pairwise within-method exact agreement of normalized relevance labels.}')
$lines.Add('\label{tab:s_pairwise_agreement}')
$lines.Add('\begin{tabular}{llr}')
$lines.Add('\toprule')
$lines.Add('Method & Run pair & Exact agreement\\')
$lines.Add('\midrule')
foreach ($row in $agreement) {
    $method = Escape-Latex $row.Method
    $pair = Escape-Latex $row.Pair
    $value = [string]::Format('{0:0.000}', $row.Value)
    $lines.Add("$method & $pair & $value\\")
}
$lines.Add('\bottomrule')
$lines.Add('\end{tabular}')
$lines.Add('\end{table}')
Set-Content -Path $pairwiseTablePath -Value $lines -Encoding UTF8

Write-Host "Generated figures in $figDir"
Write-Host "Generated supplementary tables in $suppDir"
