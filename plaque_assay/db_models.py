from typing import Any

import sqlalchemy as sql
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()  # type: Any


class NE_workflow_tracking(Base):
    __tablename__ = "NE_workflow_tracking"
    id = sql.Column(sql.Integer, primary_key=True)
    master_plate = sql.Column(sql.String(45), nullable=False)
    master_plate_id = sql.Column(sql.Integer, sql.ForeignKey("NE_master_RF_plate.id"))
    assay_plate_type = sql.Column(sql.String(45))
    start_date = sql.Column(sql.TIMESTAMP, nullable=False)
    dilution_plate_10 = sql.Column(sql.String(45))
    dilution_plate_40 = sql.Column(sql.String(45))
    dilution_plate_160 = sql.Column(sql.String(45))
    dilution_plate_640 = sql.Column(sql.String(45))
    dilutions_complete = sql.Column(sql.Integer, default=0)
    no_of_variants = sql.Column(sql.Integer, default=0)
    assays_complete = sql.Column(sql.Integer, default=0)
    sent_cl3 = sql.Column(sql.Integer, default=0)
    return_cl3 = sql.Column(sql.Integer, default=0)
    staining = sql.Column(sql.Integer, default=0)
    final_results_upload = sql.Column(sql.TIMESTAMP)
    end_date = sql.Column(sql.TIMESTAMP)
    status = sql.Column(sql.String(45))
    workflow_id = sql.Column(sql.Integer, nullable=False)


class NE_master_RF_plate(Base):
    __tablename__ = "NE_master_RF_plate"
    id = sql.Column(sql.Integer, primary_key=True)
    master_plate = sql.Column(sql.String(45), nullable=False, unique=True)
    left_rf_barcode = sql.Column(sql.String(45), nullable=False)
    right_rf_barcode = sql.Column(sql.String(45), nullable=False)
    date_uploaded = sql.Column(sql.TIMESTAMP, nullable=False)
    heat_inactivation = sql.Column(sql.Boolean, default=False)
    heat_inact_operator = sql.Column(sql.String(45))
    heat_inact_date = sql.Column(sql.TIMESTAMP, nullable=False)


class NE_raw_index(Base):
    __tablename__ = "NE_raw_index"
    id = sql.Column(sql.Integer, primary_key=True)
    row = sql.Column(sql.Integer, nullable=False)
    column = sql.Column(sql.Integer, nullable=False)
    field = sql.Column(sql.Integer, nullable=False)
    channel_id = sql.Column(sql.Integer, nullable=False)
    channel_name = sql.Column(sql.String(45), nullable=False)
    channel_type = sql.Column(sql.String(45), nullable=False)
    url = sql.Column(sql.String(240), nullable=False)
    image_resolutionx = sql.Column(sql.String(45), nullable=False)
    image_resolutiony = sql.Column(sql.String(45), nullable=False)
    image_sizex = sql.Column(sql.Integer, nullable=False)
    image_sizey = sql.Column(sql.Integer, nullable=False)
    positionx = sql.Column(sql.DECIMAL(30, 30))
    positiony = sql.Column(sql.DECIMAL(30, 30))
    time_stamp = sql.Column(sql.String(45), nullable=False)
    plate_barcode = sql.Column(sql.String(45), nullable=False)
    workflow_id = sql.Column(
        sql.Integer, sql.ForeignKey("NE_workflow_tracking.workflow_id")
    )
    variant = sql.Column(sql.String(45))


class NE_raw_results(Base):
    __tablename__ = "NE_raw_results"
    id = sql.Column(sql.Integer, primary_key=True)
    row = sql.Column(sql.Integer, nullable=False)
    column = sql.Column(sql.Integer, nullable=False)
    VPG_area_mean = sql.Column(sql.DECIMAL(30, 10))
    VPG_intensity_mean_per_well = sql.Column(sql.DECIMAL(30, 15))
    VPG_intensity_stddev_mean_per_well = sql.Column(sql.DECIMAL(30, 15))
    VPG_intensity_median_per_well = sql.Column(sql.Integer)
    VPG_intensity_sum_per_well = sql.Column(sql.BIGINT)
    cells_intensity_mean_per_well = sql.Column(sql.DECIMAL(30, 15))
    cells_intensity_stddev_mean_per_well = sql.Column(sql.DECIMAL(30, 15))
    cells_intensity_median_mean_per_well = sql.Column(sql.Integer)
    cells_intensity_sum_mean_per_well = sql.Column(sql.Integer)
    cells_image_region_area_mean_per_well = sql.Column(sql.DECIMAL(30, 15))
    normalised_plaque_area = sql.Column(sql.DECIMAL(30, 20))
    normalised_plaque_intensity = sql.Column(sql.DECIMAL(30, 20))
    number_analyzed_fields = sql.Column(sql.Integer)
    dilution = sql.Column(sql.DECIMAL(30, 20))
    well = sql.Column(sql.String(45), nullable=False)
    plate_num = sql.Column(sql.Integer)
    plate_barcode = sql.Column(sql.String(45), nullable=False)
    background_subtracted_plaque_area = sql.Column(sql.DECIMAL(30, 30))
    workflow_id = sql.Column(
        sql.Integer, sql.ForeignKey("NE_workflow_tracking.workflow_id")
    )
    variant = sql.Column(sql.String(45))


class NE_normalized_results(Base):
    __tablename__ = "NE_normalized_results"
    id = sql.Column(sql.Integer, primary_key=True)
    well = sql.Column(sql.String(45), nullable=False)
    row = sql.Column(sql.String(45), nullable=False)
    column = sql.Column(sql.String(45), nullable=False)
    dilution = sql.Column(sql.DECIMAL(30, 15))
    plate_barcode = sql.Column(sql.String(45), nullable=False)
    background_subtracted_plaque_area = sql.Column(sql.DECIMAL(30, 25))
    percentage_infected = sql.Column(sql.DECIMAL(30, 15))
    workflow_id = sql.Column(
        sql.Integer, sql.ForeignKey("NE_workflow_tracking.workflow_id")
    )
    variant = sql.Column(sql.String(45))


class NE_final_results(Base):
    __tablename__ = "NE_final_results"
    id = sql.Column(sql.Integer, primary_key=True)
    master_plate = sql.Column(sql.String(45))
    well = sql.Column(sql.String(45), nullable=False)
    ic50 = sql.Column(sql.DECIMAL(30, 15))
    status = sql.Column(sql.String(45))
    experiment = sql.Column(sql.String(45))
    workflow_id = sql.Column(
        sql.Integer, sql.ForeignKey("NE_workflow_tracking.workflow_id")
    )
    variant = sql.Column(sql.String(45))


class NE_failed_results(Base):
    __tablename__ = "NE_failed_results"
    id = sql.Column(sql.Integer, primary_key=True)
    failure_type = sql.Column(sql.String(45), nullable=False)
    plate = sql.Column(sql.String(45), nullable=False)
    well = sql.Column(sql.String(45), nullable=False)
    failure_reason = sql.Column(sql.Text, nullable=False)
    experiment = sql.Column(sql.String(45))
    workflow_id = sql.Column(
        sql.Integer, sql.ForeignKey("NE_workflow_tracking.workflow_id")
    )
    variant = sql.Column(sql.String(45))


class NE_model_parameters(Base):
    __tablename__ = "NE_model_parameters"
    id = sql.Column(sql.Integer, primary_key=True)
    well = sql.Column(sql.String(45), nullable=False)
    param_top = sql.Column(sql.DECIMAL(20, 15))
    param_bottom = sql.Column(sql.DECIMAL(20, 15))
    param_ec50 = sql.Column(sql.DECIMAL(20, 15))
    param_hillslope = sql.Column(sql.DECIMAL(20, 15))
    mean_squared_error = sql.Column(sql.DECIMAL(20, 15))
    workflow_id = sql.Column(
        sql.Integer, sql.ForeignKey("NE_workflow_tracking.workflow_id")
    )
    variant = sql.Column(sql.String(45))


class NE_assay_plate_tracker_384(Base):
    __tablename__ = "NE_assay_plate_tracker_384"
    id = sql.Column(sql.Integer, primary_key=True)
    assay_plate_384 = sql.Column(sql.String(45))
    ap_10 = sql.Column(sql.String(45))
    ap_40 = sql.Column(sql.String(45))
    ap_160 = sql.Column(sql.String(45))
    ap_640 = sql.Column(sql.String(45))
    workflow_id = sql.Column(
        sql.Integer, sql.ForeignKey("NE_workflow_tracking.workflow_id")
    )


class NE_available_strains(Base):
    __tablename__ = "NE_available_strains"
    id = sql.Column(sql.Integer, primary_key=True)
    mutant_strain = sql.Column(sql.String(45))
    plate_id_1 = sql.Column(sql.String(5))
    plate_id_2 = sql.Column(sql.String(5))


class NE_reporter_plate_status(Base):
    __tablename__ = "NE_reporter_plate_status"
    id = sql.Column(sql.Integer, primary_key=True)
    workflow_id = sql.Column(sql.Integer, nullable=False)
    variant = sql.Column(sql.String(45), nullable=False)
    status = sql.Column(sql.String(45), nullable=False)
    reporter = sql.Column(sql.String(45))
    updated_at = sql.Column(sql.TIMESTAMP)
    reason = sql.Column(sql.TEXT)


class NE_titration_workflow_tracking(Base):
    __tablename__ = "NE_titration_workflow_tracking"
    id = sql.Column(sql.Integer, primary_key=True)
    variant = sql.Column(sql.String(45))
    plate_1 = sql.Column(sql.String(45))
    plate_2 = sql.Column(sql.String(45))
    virus_batch_no = sql.Column(sql.String(45))
    start_date = sql.Column(sql.TIMESTAMP, nullable=False)
    end_date = sql.Column(sql.TIMESTAMP)
    complete_operator = sql.Column(sql.String(45))
    status = sql.Column(sql.String(45))
    workflow_id = sql.Column(sql.Integer, nullable=False)


class NE_virus_titration_normalised_results(Base):
    __tablename__ = "NE_virus_titration_normalised_results"
    id = sql.Column(sql.Integer, primary_key=True)
    plaque_area = sql.Column(sql.DECIMAL(30, 20))
    normalised_plaque_area = sql.Column(sql.DECIMAL(30, 20))
    background_subtracted_plaque_area = sql.Column(sql.DECIMAL(30, 30))
    percentage_infected = sql.Column(sql.DECIMAL(20, 15))
    dilution = sql.Column(sql.Integer, nullable=False)
    well = sql.Column(sql.String(3), nullable=False)
    plate_barcode = sql.Column(sql.String(9))
    workflow_id = sql.Column(sql.Integer, nullable=False)


class NE_virus_titration_model_parameters(Base):
    __tablename__ = "NE_virus_titration_model_parameters"
    id = sql.Column(sql.Integer, primary_key=True)
    dilution = sql.Column(sql.Integer, nullable=False)
    param_top = sql.Column(sql.DECIMAL(20, 15))
    param_bottom = sql.Column(sql.DECIMAL(20, 15))
    param_ec50 = sql.Column(sql.DECIMAL(20, 15))
    param_hillslope = sql.Column(sql.DECIMAL(20, 15))
    mean_squared_error = sql.Column(sql.DECIMAL(20, 15))
    workflow_id = sql.Column(sql.Integer, nullable=False)


class NE_virus_titration_final_results(Base):
    __tablename__ = "NE_virus_titration_final_results"
    id = sql.Column(sql.Integer, primary_key=True)
    dilution = sql.Column(sql.Integer, nullable=False)
    ic50 = sql.Column(sql.DECIMAL(30, 15))
    status = sql.Column(sql.String(45))
    workflow_id = sql.Column(sql.Integer, nullable=False)
