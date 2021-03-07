
# require 'rubygems'
# require 'sitemap_generator'

SitemapGenerator::Sitemap.default_host = 'https://primer-monitor.neb.com'
SitemapGenerator::Sitemap.create do
  add '/', changefreq: 'daily'
  add organisms_path, changefreq: 'yearly'
  add lineages_path, changefreq: 'weekly'
  add about_path, changefreq: 'yearly'
  add primer_sets_path, changefreq: 'daily'
  PrimerSet.find_each do |ps|
    add primer_set_path(ps), lastmod: ps.updated_at
  end
end
SitemapGenerator::Sitemap.ping_search_engines # Not needed if you use the rake tasks
