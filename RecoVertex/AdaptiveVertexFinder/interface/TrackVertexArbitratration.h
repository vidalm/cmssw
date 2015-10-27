#ifndef TrackVertexArbitration_H
#define TrackVertexArbitration_H
#include <memory>
#include <set>
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/CandidatePtrTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
//#include "DataFormats/Math/interface/deltaR.h"
#define VTXDEBUG
template <class VTX>
class TrackVertexArbitration{
public:
TrackVertexArbitration(const edm::ParameterSet &params);
std::vector<VTX> trackVertexArbitrator(
edm::Handle<reco::BeamSpot> &beamSpot,
const reco::Vertex &pv,
std::vector<reco::TransientTrack> & selectedTracks,
std::vector<VTX> & secondaryVertices
);
private:
bool trackFilterArbitrator(const reco::TransientTrack &track) const;
edm::InputTag primaryVertexCollection;
edm::InputTag secondaryVertexCollection;
edm::InputTag trackCollection;
edm::InputTag beamSpotCollection;
double dRCut;
double	distCut;
double	sigCut;
double	dLenFraction;
double fitterSigmacut; // = 3.0;
double fitterTini; // = 256.;
double fitterRatio;// = 0.25;
int trackMinLayers;
double	trackMinPt;
int	trackMinPixels;
};
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
template <class VTX>
TrackVertexArbitration<VTX>::TrackVertexArbitration(const edm::ParameterSet &params) :
primaryVertexCollection (params.getParameter<edm::InputTag>("primaryVertices")),
  secondaryVertexCollection (params.getParameter<edm::InputTag>("secondaryVertices")),
  trackCollection (params.getParameter<edm::InputTag>("tracks")),
  beamSpotCollection (params.getParameter<edm::InputTag>("beamSpot")),
  dRCut (params.getParameter<double>("dRCut")),
  distCut (params.getParameter<double>("distCut")),
  sigCut (params.getParameter<double>("sigCut")),
  dLenFraction (params.getParameter<double>("dLenFraction")),
  fitterSigmacut (params.getParameter<double>("fitterSigmacut")),
  fitterTini (params.getParameter<double>("fitterTini")),
  fitterRatio (params.getParameter<double>("fitterRatio")),
  trackMinLayers (params.getParameter<int32_t>("trackMinLayers")),
  trackMinPt (params.getParameter<double>("trackMinPt")),
  trackMinPixels (params.getParameter<int32_t>("trackMinPixels"))
{
  dRCut*=dRCut;
}
template <class VTX>
bool TrackVertexArbitration<VTX>::trackFilterArbitrator(const reco::TransientTrack &track) const
{
  if(!track.isValid())
    return false;
  if (track.track().hitPattern().trackerLayersWithMeasurement() < trackMinLayers)
    return false;
  if (track.track().pt() < trackMinPt)
    return false;
  if (track.track().hitPattern().numberOfValidPixelHits() < trackMinPixels)
    return false;
  return true;
}
float trackWeight(const reco::Vertex & sv, const reco::TransientTrack &track)
{
return sv.trackWeight(track.trackBaseRef());
}
float trackWeight(const reco::VertexCompositePtrCandidate & sv, const reco::TransientTrack &tt)
{
  const reco::CandidatePtrTransientTrack* cptt = dynamic_cast<const reco::CandidatePtrTransientTrack*>(tt.basicTransientTrack());
  if ( cptt==0 )
    edm::LogError("DynamicCastingFailed") << "Casting of TransientTrack to CandidatePtrTransientTrack failed!";
  else
    {
      const	reco::CandidatePtr & c=cptt->candidate();
      if(std::find(sv.daughterPtrVector().begin(),sv.daughterPtrVector().end(),c)!= sv.daughterPtrVector().end())
	return 1.0;
      else
	return 0;
    }
  return 0;
}
template <class VTX>
std::vector<VTX> TrackVertexArbitration<VTX>::trackVertexArbitrator(
edm::Handle<reco::BeamSpot> &beamSpot,
const reco::Vertex &pv,
std::vector<reco::TransientTrack> & selectedTracks,
std::vector<VTX> & secondaryVertices)
{
  using namespace reco;
  //std::cout << "PV: " << pv.position() << std::endl;
  VertexDistance3D dist;
  VertexDistance3D mydist;
  AdaptiveVertexFitter theAdaptiveFitter(
					 GeometricAnnealing(fitterSigmacut, fitterTini, fitterRatio),
					 DefaultLinearizationPointFinder(),
					 KalmanVertexUpdator<5>(),
					 KalmanVertexTrackCompatibilityEstimator<5>(),
					 KalmanVertexSmoother() );
  std::vector<VTX> recoVertices;
  VertexDistance3D vdist;
  GlobalPoint ppv(pv.position().x(),pv.position().y(),pv.position().z());
  std::map<unsigned int, Measurement1D> cachedIP;
  for(typename std::vector<VTX>::const_iterator sv = secondaryVertices.begin();
      sv != secondaryVertices.end(); ++sv) {
    GlobalPoint ssv(sv->position().x(),sv->position().y(),sv->position().z());
    GlobalVector flightDir = ssv-ppv;
    // std::cout << "Vertex : " << sv-secondaryVertices->begin() << " " << sv->position() << std::endl;
    Measurement1D dlen= vdist.distance(pv,VertexState(RecoVertex::convertPos(sv->position()),RecoVertex::convertError(sv->error())));
    std::vector<reco::TransientTrack> selTracks;
    for(unsigned int itrack = 0; itrack < selectedTracks.size(); itrack++){
      TransientTrack & tt = (selectedTracks)[itrack];
      if (!trackFilterArbitrator(tt)) continue;
      tt.setBeamSpot(*beamSpot);
      float w = trackWeight(*sv,tt);
      Measurement1D ipv;
      if(cachedIP.find(itrack)!=cachedIP.end()) {
	ipv=cachedIP[itrack];
      } else {
	std::pair<bool,Measurement1D> ipvp = IPTools::absoluteImpactParameter3D(tt,pv);
	cachedIP[itrack]=ipvp.second;
	ipv=ipvp.second;
      }
      AnalyticalImpactPointExtrapolator extrapolator(tt.field());
      TrajectoryStateOnSurface tsos = extrapolator.extrapolate(tt.impactPointState(), RecoVertex::convertPos(sv->position()));
      if(! tsos.isValid()) continue;
      GlobalPoint refPoint = tsos.globalPosition();
      GlobalError refPointErr = tsos.cartesianError().position();
      Measurement1D isv = dist.distance(VertexState(RecoVertex::convertPos(sv->position()),RecoVertex::convertError(sv->error())),VertexState(refPoint, refPointErr));
      float dR = Geom::deltaR2(flightDir,tt.track()); //.eta(), flightDir.phi(), tt.track().eta(), tt.track().phi());
      if( w > 0 || ( isv.significance() < sigCut && isv.value() < distCut && isv.value() < dlen.value()*dLenFraction ) )
	{
	  if(( isv.value() < ipv.value() ) && isv.value() < distCut && isv.value() < dlen.value()*dLenFraction
	     && dR < dRCut )
	    {

	      //selTracks.push_back(tt);
	      // Check if the track distance to sv is smaller than to another secondary vertices
	      Measurement1D minDist=isv;
	      GlobalPoint ssvmin=ssv;

	      for(typename std::vector<VTX>::const_iterator svj = secondaryVertices.begin();
		  svj != secondaryVertices.end(); ++svj) {
		if (sv == svj) continue;
		Measurement1D svdist = mydist.distance(VertexState(RecoVertex::convertPos(sv->position()),RecoVertex::convertError(sv->error())),VertexState(RecoVertex::convertPos(svj->position()),RecoVertex::convertError(svj->error())));
		if(svdist.significance() > 50){
		  TrajectoryStateOnSurface tsosj = extrapolator.extrapolate(tt.impactPointState(), RecoVertex::convertPos(svj->position()));
		  if(! tsosj.isValid()) continue;
		  GlobalPoint refPointj = tsosj.globalPosition();
		  GlobalError refPointErrj = tsosj.cartesianError().position();
		  Measurement1D isvj = dist.distance(VertexState(RecoVertex::convertPos(svj->position()),RecoVertex::convertError(svj->error())),VertexState(refPointj, refPointErrj));
		  if (isvj.significance() < minDist.significance()){
		    minDist=isvj;
		    GlobalPoint ssvj(svj->position().x(),svj->position().y(),svj->position().z());
		    ssvmin=ssvj;
		  }
		}
	      }
	      GlobalVector flightDirmin = ssvmin-ppv;
	      float dRmin = Geom::deltaR2(flightDirmin,tt.track());
	      if (minDist.significance() >= 0.01*isv.significance() || dRmin > dR) selTracks.push_back(tt);
	    
	      
#ifdef VTXDEBUG
	      if (!(minDist.significance() >= 0.01*isv.significance() || dRmin > dR) ){ 
		std::cout <<"NOT Pushed back"<< std::endl;
		std::cout <<"Vertex orig "<<sv-secondaryVertices.begin() <<  std::endl;
		std::cout <<"svip:     " << isv.significance() << " " << isv.value()<<" "<<dR <<std::endl;
		std::cout <<"min_svip: " << minDist.significance() << " " << minDist.value()<<" " <<dRmin <<std::endl;
		
	      }
#endif
	      
	    } else
	    {
	      //add also the tracks used in previous fitting that are still closer to Sv than Pv
	      if(w > 0.5 && isv.value() <= ipv.value() && dR < dRCut) {
		
		// Check if the track distance to sv is smaller than to another secondary vertices
		Measurement1D minDist=isv;
		GlobalPoint ssvmin=ssv;
		
		for(typename std::vector<VTX>::const_iterator svj = secondaryVertices.begin();
		    svj != secondaryVertices.end(); ++svj) {
		  if (sv == svj) continue;
		  Measurement1D svdist = mydist.distance(VertexState(RecoVertex::convertPos(sv->position()),RecoVertex::convertError(sv->error())),VertexState(RecoVertex::convertPos(svj->position()),RecoVertex::convertError(svj->error())));
		  
		  if(svdist.significance() > 50){
		    TrajectoryStateOnSurface tsosj = extrapolator.extrapolate(tt.impactPointState(), RecoVertex::convertPos(svj->position()));
		    if(! tsosj.isValid()) continue;
		    GlobalPoint refPointj = tsosj.globalPosition();
		    GlobalError refPointErrj = tsosj.cartesianError().position();
		    Measurement1D isvj = dist.distance(VertexState(RecoVertex::convertPos(svj->position()),RecoVertex::convertError(svj->error())),VertexState(refPointj, refPointErrj));
		    if (isvj.significance() < minDist.significance()){
		      minDist=isvj;
		      GlobalPoint ssvj(svj->position().x(),svj->position().y(),svj->position().z());
		      ssvmin=ssvj;
		    }
		  }
		}
		GlobalVector flightDirmin = ssvmin-ppv;
		float dRmin = Geom::deltaR2(flightDirmin,tt.track());
		if (minDist.significance() >= 0.01*isv.significance() || dRmin > dR) selTracks.push_back(tt);
		// selTracks.push_back(tt);
#ifdef VTXDEBUG
		if (! (minDist.significance() >= 0.01*isv.significance()|| dRmin > dR) ){ 
		  std::cout <<"NOT Pushed back 2ND STEP"<< std::endl;
		  std::cout <<"Vertex orig "<<sv-secondaryVertices.begin() <<  std::endl;
		  std::cout <<"svip:     " << isv.significance() << " " << isv.value()<< std::endl;
		  std::cout <<"min_svip: " << minDist.significance() << " " << minDist.value()<<" " <<dRmin <<std::endl;
		}
#endif

	      }
	    }
	}
    }
    if(selTracks.size() >= 2)
      {
	TransientVertex singleFitVertex;
	singleFitVertex = theAdaptiveFitter.vertex(selTracks,ssv);
	if(singleFitVertex.isValid()) { recoVertices.push_back(VTX(singleFitVertex));
#ifdef VTXDEBUG	  
	if(!singleFitVertex.isValid()) std::cout <<"Vertex fit NOT VALID"<< std::endl;
#endif	  
	}
      }
  }
  return recoVertices;
}
#endif
