import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import path from 'path'

const mathjax = mathjaxPlugin()

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  ignoreDeadLinks: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  markdown: {
    config(md) {
      md.use(tabsMarkdownPlugin);
      md.use(footnote);
      mathjax.markdownConfig(md);
    },
    codeTransformers: [juliaReplTransformer()],
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  vite: {
    plugins: [
      mathjax.vitePlugin,
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  themeConfig: {
    outline: 'deep',
    siteTitle: false,
    logo: {
      light: '/logo.svg',
      dark: '/logo-dark.svg'
    },
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: {
      '/tutorial/': [
        {
          text: 'Tutorials',
          collapsed: false,
          items: [
            { text: 'Free Induction Decay', link: '/tutorial/01-FID' },
            { text: 'Small Tip Angle Approximation', link: '/tutorial/02-SmallTipApproximation' },
            { text: 'Chemical Shift in an EPI sequence', link: '/tutorial/03-ChemicalShiftEPI' },
            { text: 'Slice-Selective Acquisition of 3D Phantom', link: '/tutorial/04-3DSliceSelective' },
            { text: "Patient's Motion During Acquisition", link: '/tutorial/05-SimpleMotion' },
            { text: 'Diffusion-induced Signal Attenuation', link: '/tutorial/06-DiffusionMotion' },
            { text: 'Cardiac Cine MRI with Arrhythmias', link: '/tutorial/07-RRVariability' },
            { text: 'Using Labels to reconstruct multi-slice / multi-contrast sequences', link: '/tutorial/07-label' },
          ]
        },
        {
          text: 'Reproducible Tutorials',
          collapsed: false,
          items: [
            { text: 'Understanding basic MRI sequences', link: '/tutorial-pluto/01-gradient-echo-spin-echo' },
            { text: 'Low-Field CMRA Optimization', link: '/tutorial-pluto/02-low-field-cmra-optimization' },
            { text: 'Low-Field BOOST Optimization', link: '/tutorial-pluto/03-low-field-boost-optimization' },
          ]
        }
      ],
      '/tutorial-pluto/': [
        {
          text: 'Tutorials',
          collapsed: false,
          items: [
            { text: 'Free Induction Decay', link: '/tutorial/01-FID' },
            { text: 'Small Tip Angle Approximation', link: '/tutorial/02-SmallTipApproximation' },
            { text: 'Chemical Shift in an EPI sequence', link: '/tutorial/03-ChemicalShiftEPI' },
            { text: 'Slice-Selective Acquisition of 3D Phantom', link: '/tutorial/04-3DSliceSelective' },
            { text: "Patient's Motion During Acquisition", link: '/tutorial/05-SimpleMotion' },
            { text: 'Diffusion-induced Signal Attenuation', link: '/tutorial/06-DiffusionMotion' },
            { text: 'Cardiac Cine MRI with Arrhythmias', link: '/tutorial/07-RRVariability' },
            { text: 'Using Labels to reconstruct multi-slice / multi-contrast sequences', link: '/tutorial/07-label' },
          ]
        },
        {
          text: 'Reproducible Tutorials',
          collapsed: false,
          items: [
            { text: 'Understanding basic MRI sequences', link: '/tutorial-pluto/01-gradient-echo-spin-echo' },
            { text: 'Low-Field CMRA Optimization', link: '/tutorial-pluto/02-low-field-cmra-optimization' },
            { text: 'Low-Field BOOST Optimization', link: '/tutorial-pluto/03-low-field-boost-optimization' },
          ]
        }
      ],
      '/how-to/': [
        {
          text: 'How To',
          collapsed: false,
          items: [
            { text: 'Getting Started', link: '/how-to/1-getting-started' },
            { text: 'Use Koma UI', link: '/how-to/2-1-use-koma-ui' },
            { text: 'Use Koma Notebooks', link: '/how-to/2-2-use-koma-notebooks' },
            { text: 'Use Koma Scripts', link: '/how-to/2-3-use-koma-scripts' },
            { text: 'Create Your Own Phantom', link: '/how-to/3-create-your-own-phantom' },
            { text: 'Create Your Own Sequence', link: '/how-to/3-create-your-own-sequence' },
            { text: 'Run Distributed Simulations', link: '/how-to/4-run-distributed-simulations' },
            { text: 'Contribute to Koma', link: '/how-to/5-contribute-to-koma' },
          ]
        }
      ],
      '/explanation/': [
        {
          text: 'Explanations',
          collapsed: false,
          items: [
            { text: 'Phantom', link: '/explanation/1-phantom' },
            { text: 'Motion', link: '/explanation/2-motion' },
            { text: 'Phantom Format', link: '/explanation/3-phantom-format' },
            { text: 'Sequence', link: '/explanation/4-sequence' },
            { text: 'Sequence Events', link: '/explanation/5-seq-events' },
            { text: 'Simulation', link: '/explanation/6-simulation' },
            { text: 'GPU Explanation', link: '/explanation/7-gpu-explanation' },
          ]
        }
      ],
      '/reference/': [
        {
          text: 'Reference',
          collapsed: false,
          items: [
            { text: 'API', link: '/reference/1-api' },
            { text: 'KomaMRIBase', link: '/reference/2-koma-base' },
            { text: 'KomaMRICore', link: '/reference/3-koma-core' },
            { text: 'KomaMRIFiles', link: '/reference/4-koma-files' },
            { text: 'KomaMRIPlots', link: '/reference/5-koma-plots' },
            { text: 'KomaMRI', link: '/reference/6-koma-mri' },
          ]
        }
      ],
    },
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: "slack", link: "https://julialang.org/slack/" },
    ],
    footer: {
      message:
        'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a>, and <a href="https://luxdl.github.io/DocumenterVitepress.jl/stable" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>Released under the MIT License. Powered by the <a href="https://www.julialang.org">Julia Programming Language</a>.<br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()} Koma Development Team.`,
    }
  }
})